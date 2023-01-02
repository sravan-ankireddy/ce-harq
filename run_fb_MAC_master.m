% Fixing the range for reproducibility 
rng("default");

if ~exist('global_settings','var')
    % Code parameters
    targetCodeRate = 0.5;
    
    N = 400;
    code_type = "Conv";
    K = round(N*targetCodeRate);
    R = targetCodeRate;
    combining_scheme = "CC";
    dec_type = "hard";

    SNRdB_low = -10;
    SNRdB_high = -5;
    SNRdB_step = 0.5;
    SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
    
    num_SNRdB = length(SNRdB_vec);
    err_thr = 0.03;
    err_thr_ada_list = err_thr*ones(size(SNRdB_vec));
    err_thr_ada_list_est = 0;
    max_rounds = 10;
    modulation = 'BPSK';
    k = bits_per_symbol(modulation);
    M = 2^k;
    mod_approx = 0;
    comm_mod = 1;
    
    err_thr_ada_scheme = "opt";
    channel = "awgn";
    
    % num frames
    nOut = 10;
    nMiniFrames = 1000;
    nMinFerr = 100;

    res_folder_all = "results/all";

    process_data_fb = 0;
end

% Error Stats
BER_vec_FB = zeros(1,num_SNRdB);
BLER_vec_FB = zeros(1,num_SNRdB);

BER_vec_pr_FB = zeros(max_rounds,num_SNRdB);
BLER_vec_pr_FB = zeros(max_rounds,num_SNRdB);

Avg_rounds_FB = zeros(1,num_SNRdB);


tic;
for i_s = 1:length(SNRdB_vec)
    SNRdB = round(SNRdB_vec(i_s),4);

    noiseVar = 1./(10.^(SNRdB/10));

    err_thr = err_thr_ada_list(i_s);

    if (err_thr_ada_scheme == "est")
        err_thr = err_thr_ada_list_est((max_rounds - 1),i_s);
    end

    % Compressor set up
    p = round(100*err_thr);
    counts = [100-p p];

    % Print status
    fprintf('Status %0.2f %% done \n', round(i_s/num_SNRdB*100));
    
    BER_FB = 0; BLER_FB = 0;   
    BER_FB_per_round = zeros(max_rounds,1); BLER_FB_per_round = zeros(max_rounds,1);
    % Count the average rounds for successful decoding
    num_ar_fb = 0;

    for i_on = 1:nOut
    
        parfor i_n = 1:nMiniFrames
            
            seed = i_n + (i_on-1)*nMiniFrames;
            
            rng(2*seed);

            % Update the counter 
            num_ar_fb = num_ar_fb + 1; 

            % Generate random message
            data = randi([0 1], K, 1);
        
            % Encoding and Rate matching
            if (code_type == "Conv")
                [dataIn, rr_len] = conv_enc(data, targetCodeRate);
            else
                bgn = bgn_select(K,R);
                dataIn = nrldpc_enc(data, R, modulation, rv, bgn, nlayers);
            end
            
            % Symbol Modulation
            if (modulation == "BPSK")
                bpskModulator = comm.BPSKModulator;
                bpskModulator.PhaseOffset = pi/4;
                txSig = bpskModulator(dataIn);
            else
                txSig = qammod(dataIn,M,'InputType','bit','UnitAveragePower',true);
            end

            % Pass through channel
            rxSig = txSig;
            if (channel == "awgn")
                rxSig = awgn(txSig,SNRdB,'measured');
            elseif (channel == "rayleigh")
                h = sqrt(rand(1)^2 + rand(1)^2);
                rxSig = awgn(h*txSig,SNRdB,'measured');
            end

            % Symbol demod
            if (modulation == "BPSK")
                if (dec_type == "hard")
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Hard decision';
                    rxLLR = bpskDemodulator(rxSig);
                else
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Approximate log-likelihood ratio';
                    rxLLR = bpskDemodulator(rxSig);
                end
            else
                if (dec_type == "hard")
                    rxLLR = qamdemod(rxSig, M, OutputType='bit', UnitAveragePower=true);
                else
                    rxLLR = qamdemod(rxSig, M, OutputType='approxllr', UnitAveragePower=true, NoiseVariance=noiseVar);
                end
            end

            % Rate recovery and Decoding
            if (code_type == "Conv")
                data_est = conv_dec(rxLLR, targetCodeRate, dec_type);
            else
                rxLLR_rr = nrRateRecoverLDPC(rxLLR, K, R, rv, modulation, nlayers, ncb, Nref);
                bgn = bgn_select(K,R);
                [data_est, crc_chk] = nrldpc_dec(rxLLR_rr, K, max_iter, bgn);
            end

            % Check for error stats
            num_err = sum(data ~= double(data_est));
            num_err_FB = num_err;
            num_err_FB_per_round = zeros(max_rounds,1);
            num_err_FB_per_round(1) = num_err;
            
            % Start retransmission if 1st round failed
            if (num_err > 0 && max_rounds > 1)
                out = retransmit_func_FB_MAC_master(channel,SNRdB,modulation,N,K,targetCodeRate,code_type,ncb,Nref,max_iter,nlayers,dec_type,data,rxLLR,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,comm_mod,mod_approx,seed);
                num_ar_fb = num_ar_fb + out.Avg_rounds_FB;
                num_err_FB = out.num_err_FB;
                num_err_FB_per_round = out.num_err_vec;
            end

            BER_FB = BER_FB + num_err_FB;

            BLER_FB = BLER_FB + (num_err_FB > 0);

            BER_FB_per_round = BER_FB_per_round + num_err_FB_per_round;

            BLER_FB_per_round = BLER_FB_per_round + (num_err_FB_per_round > 0);
        end
        if (BLER_FB > nMinFerr)
            break;
        end
    end
    % AR update
    Avg_rounds_FB(i_s) = num_ar_fb/(i_on*nMiniFrames);

    % Error Stats
    BER_vec_FB(:,i_s) = BER_FB/(K*(i_on*nMiniFrames));
    BLER_vec_FB(:,i_s) = BLER_FB/(i_on*nMiniFrames);

    BER_vec_pr_FB(:,i_s) = BER_FB_per_round/(K*(i_on*nMiniFrames));
    BLER_vec_pr_FB(:,i_s) = BLER_FB_per_round/(i_on*nMiniFrames);
end


toc;

%% 
if (err_thr_ada_scheme == "opt")
    disp("err_thr_used : opt: ");
    disp(err_thr_ada_list);
else
    disp("err_thr_used : est: ");
    disp(err_thr_ada_list_est);
end


% plotting
if (process_data_fb == 1)
    
    % BER
    % set the figure properties
    figure('Renderer','painters','Position',[1000 400 800 500]);
    
    f = semilogy(SNRdB_vec,BER_vec_pr_FB(end,:),'r-d');
    
    fs = 12;
    xlabel('SNR','FontSize',fs);
    ylabel('BER','FontSize',fs);
    grid on;
    
    codeRate = R;
    leg_FB = sprintf('FB-%s BER Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);
    
    legend(leg_FB, 'Location','southwest','FontSize',fs);
    
    title_name = sprintf('FB-%s : BER %s %d, mod. %s, Rate %.3f, errthr %.3f, max. rounds %d',code_type, combining_scheme, N, modulation, R, err_thr, max_rounds);
    title(title_name,'FontSize',fs);
    
    common_str = res_folder_all + sprintf('/BER_%s_%d_rate_%.3f_err_thr_%.3f_max_rounds_%d', code_type, N, R, err_thr, max_rounds);
    
    fig_name = common_str + ".fig";
    savefig(fig_name);
    png_name = common_str + ".png";
    saveas(f,png_name);

    % BLER
    % set the figure properties
    figure('Renderer','painters','Position',[1000 400 800 500]);
    
    f = semilogy(SNRdB_vec,BLER_vec_pr_FB(end,:),'r-d');
    
    fs = 12;
    xlabel('SNR','FontSize',fs);
    ylabel('BLER','FontSize',fs);
    grid on;
    
    codeRate = R;
    leg_FB = sprintf('FB-%s BLER Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);
    
    legend(leg_FB, 'Location','southwest','FontSize',fs);
    
    title_name = sprintf('FB-%s : BLER %s %d, mod. %s, Rate %.3f, errthr %.3f, max. rounds %d',code_type, combining_scheme, N, modulation, R, err_thr, max_rounds);
    title(title_name,'FontSize',fs);
    
    common_str = res_folder_all + sprintf('/BLER_%s_%d_rate_%.3f_err_thr_%.3f_max_rounds_%d',code_type, N, R, err_thr, max_rounds);
    
    fig_name = common_str + ".fig";
    savefig(fig_name);
    png_name = common_str + ".png";
    saveas(f,png_name);
    
    %% Avg rounds

    % set the figure properties
    figure('Renderer','painters','Position',[1000 400 800 500]);
    
    f = semilogy(SNRdB_vec,Avg_rounds_FB,'r-d');
    
    fs = 12;
    xlabel('SNR','FontSize',fs);
    ylabel('Avg. rounds','FontSize',fs);
    grid on;
    
    leg_FB = sprintf('FB-%s AR Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);
    legend(leg_FB, 'Location','southwest','FontSize',fs);
    
    title_name = sprintf('FB-%s : AR %s %d, mod. %s, Rate %.3f, errthr %.3f, max. rounds %d', code_type, combining_scheme, N, modulation, R, err_thr, max_rounds);
    title(title_name,'FontSize',fs);
    
    common_str = res_folder_all + sprintf('/AR_%s_%d_rate_%.3f_err_thr_%.3f_max_rounds_%d', code_type, N, R, err_thr, max_rounds);
    
    fig_name = common_str + ".fig";
    savefig(fig_name);
    png_name = common_str + ".png";
    saveas(f,png_name);

    %% Spectral Efficiency

    % set the figure properties
    figure('Renderer','painters','Position',[1000 400 800 500]);
    
    SE = (K*k/N) * (1 - BLER_vec_pr_FB(end,:)./Avg_rounds_FB);

    f = semilogy(SNRdB_vec,SE,'r-d');
    
    fs = 12;
    xlabel('SNR','FontSize',fs);
    ylabel('Spectral Efficiency [bits/s/Hz]','FontSize',fs);
    grid on;
    
    leg_FB = sprintf('FB-%s SE Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);
    legend(leg_FB, 'Location','southwest','FontSize',fs);
    
    title_name = sprintf('FB-%s : SE %s %d, mod. %s, Rate %.3f, errthr %.3f, max. rounds %d', code_type, combining_scheme, N, modulation, R, err_thr, max_rounds);
    title(title_name,'FontSize',fs);
    
    common_str = res_folder_all + sprintf('/SE_%s_%d_rate_%.3f_err_thr_%.3f_max_rounds_%d', code_type, N, R, err_thr, max_rounds);
    
    fig_name = common_str + ".fig";
    savefig(fig_name);
    png_name = common_str + ".png";
    saveas(f,png_name);
    
    ber_data = BER_vec_pr_FB;
    bler_data = BLER_vec_pr_FB;
    ar_data = Avg_rounds_FB;
    snr_data = SNRdB_vec;
    
    % Store data 
    data_file_name = res_folder_all + sprintf('/fb_data_%s_%d_rate_%.3f_err_thr_%.3f_max_rounds_%d.mat', code_type, N, R, err_thr, max_rounds);
    save(data_file_name,'ber_data','bler_data','ar_data','snr_data','err_thr');
end