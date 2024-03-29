% Fixing the range for reproducibility 
rng("default");

% Simulation settings : Defaults
if ~exist('global_settings','var')
    % Code parameters
    targetCodeRate = 0.8;
    nPRB = 20; % Vary this to change the code length
    nlayers = 1;
    NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice
    N = nPRB*NREPerPRB;
    max_iter = 6; % default is 8 in MATLAB
    min_ber = 1e-5;
    SNRdB_low = -8;
    SNRdB_high = 10;
    SNRdB_step = 0.2;
    SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
    num_SNRdB = length(SNRdB_vec);
    nFrames = 10e1;
    err_thr = 0.05;
    err_thr_ada_list = err_thr*ones(size(SNRdB_vec));
    max_rounds = 10;
    modulation = 'QPSK';
    mod_approx = 0;
    weighted_mean = 0;
    qam_mod = 0;
    min_frame_error = 0;
    res_folder_all = "bler_data/fb";
    process_data_fb = 0;
    err_thr_ada_scheme = "opt";
    generate_harq_lut_data == 0;
end

if (generate_harq_lut_data == 1)
    err_thr = 0;
    
    % use nrTBS to get K,R
    tbs = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetCodeRate); %no. bits in transportBlock

    % CodeLen based on PRB settings and modulation
    N = nPRB*NREPerPRB*M;
    K = tbs;
    R = K/N;
end

% base graph selection based on rate and K
bgn = bgn_select(K,R);

rv = rvSeq(1);

% Error Stats
BER_vec_FB = zeros(1,num_SNRdB);
BLER_vec_FB = zeros(1,num_SNRdB);

BER_vec_pr_FB = zeros(max_rounds,num_SNRdB);
BLER_vec_pr_FB = zeros(max_rounds,num_SNRdB);

% After nth round, store the probability of error corresponding to sparsity of error in round n-1
err_stats_vec_pr_FB = zeros(max_rounds-1,round(err_thr*N),num_SNRdB);

% After nth round, store the count of count of errors with spartsity s
err_sparsity_vec_pr_FB = zeros(max_rounds,N,num_SNRdB);

Avg_rounds_FB = zeros(size(SNRdB_vec));

% Retransmission params
fb_scheme = "HARQ";
harq_for_fb = 1;

opts.MappingType = 'signed';
opts.OutputDataType = 'double';

tic;
for i_s = 1:length(SNRdB_vec)
    SNRdB = round(SNRdB_vec(i_s),4);

    noiseVar = 1./(10.^(SNRdB/10));

    err_thr = err_thr_ada_list(i_s);

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
            bgn = bgn_select(K,R);
            dataIn = nrldpc_enc(data, R, modulation, rv, bgn, nlayers);
            
            % Scrambling, TS 38.211 Section 7.3.1.1
            nid = 1; rnti = 1; cwi = 1;
            sc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(dataIn));
            dataIn_sc = xor(dataIn,sc_seq);

            % Symbol Modulation
            if (qam_mod == 1)
                % Reshape data into binary k-tuples, k = log2(M), and convert to
                % integers
                txData = reshape(dataIn_sc,length(dataIn)/M,M);
                txDataSym = bi2de(txData);
                txSig = qammod(txDataSym,2^M,'bin','UnitAveragePower',1);
            else
                txSig = nrSymbolModulate(dataIn_sc,modulation);
            end
                
            % Pass through AWGN channel
            rxSig = awgn(txSig,SNRdB);

            % Symbol demod
            if (qam_mod == 1)
                if (mod_approx == 0)
                    rxLLR = qamdemod(rxSig,2^M,'OutputType','approxllr');
                else
                    rxLLR = 1 - 2*double(qamdemod(rxSig,2^M,'OutputType','bit','UnitAveragePower',1));
                end
            else
                if (mod_approx == 0)
                    noiseVar = 1./(10.^(SNRdB/10));
                    rxLLR = nrSymbolDemodulate(rxSig,modulation,noiseVar);
                else
                    rxLLR = 1 - 2*double(nrSymbolDemodulate(rxSig,modulation,'DecisionType','hard'));
                end
            end

            % Descrambling, inverse of TS 38.211 Section 7.3.1.1
            dsc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(rxLLR),opts);
            rxLLR_dsc = rxLLR .* dsc_seq;

            % Rate recovery
            rxLLR_dsc_rr = nrRateRecoverLDPC(rxLLR_dsc, K, R, rv, modulation, nlayers, ncb, Nref);
            % Decoding
            bgn = bgn_select(K,R);
            [data_est, crc_chk] = nrldpc_dec(rxLLR_dsc_rr, K, max_iter, bgn);

            % Check for error stats
            num_err = sum(data ~= double(data_est));
            num_err_FB = num_err;
            num_err_FB_per_round = zeros(max_rounds,1);
            num_err_FB_per_round(1) = num_err;
            
            % Start retransmission if 1st round failed
            if (crc_chk > 0 && max_rounds > 1)
                out = retransmit_func_FB(SNRdB,modulation,max_iter,rvSeq,nlayers,nPRB,NREPerPRB,N,K,R,data,rxLLR_dsc_rr,ncb,Nref,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,qam_mod,mod_approx,seed);
                num_ar_fb = num_ar_fb + out.Avg_rounds_FB;
                num_err_FB = out.num_err_FB;
                num_err_FB_per_round = out.num_err_vec;
            end

            BER_FB = BER_FB + num_err_FB;

            BLER_FB = BLER_FB + (num_err_FB > 0);

            BER_FB_per_round = BER_FB_per_round + num_err_FB_per_round;

            BLER_FB_per_round = BLER_FB_per_round + (num_err_FB_per_round > 0);
        end

        if (BLER_FB > min_err && i_on*nMiniFrames > min_blocks)
            break;
        end
    end
    % AR update
    Avg_rounds_FB(i_s) = Avg_rounds_FB(i_s) + num_ar_fb;

    % Error Stats
    BER_vec_FB(:,i_s) = BER_FB/(K*(i_on*nMiniFrames));
    BLER_vec_FB(:,i_s) = BLER_FB/(i_on*nMiniFrames);

    BER_vec_pr_FB(:,i_s) = BER_FB_per_round/(K*(i_on*nMiniFrames));
    BLER_vec_pr_FB(:,i_s) = BLER_FB_per_round/(i_on*nMiniFrames);
end

% Avg rounds per transmission
Avg_rounds_FB = Avg_rounds_FB/nFrames;

toc;

if (process_data_harq == 1)
    BER_vec_pr_HARQ = BER_vec_pr_FB;
    BLER_vec_pr_HARQ = BLER_vec_pr_FB;
    Avg_rounds_HARQ = Avg_rounds_FB;

	actualCodeRate_list(i_tr) = R;

	% current all rounds
	err_data_all_rates(1,i_tr,:,:) = BER_vec_pr_HARQ;
	err_data_all_rates(2,i_tr,:,:) = BLER_vec_pr_HARQ;

	ar_data_all_rates(i_tr,:) = Avg_rounds_HARQ;
	
	% Only for current rate 
	err_data = squeeze(err_data_all_rates(:,i_tr,:,:));
	ar_data = squeeze(ar_data_all_rates(i_tr,:));

	 % set the figure properties
	figure('Renderer','painters','Position',[1000 400 800 500]);
		
	f = semilogy(SNRdB_vec,BER_vec_pr_HARQ(end,:),'b-o');
	
	fs = 12;
	xlabel('SNR','FontSize',fs);
	ylabel('BER','FontSize',fs);
	
	codeRate = R;
	leg_HARQ = sprintf('HARQ-%s BER Rate %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
	
	legend(leg_HARQ, 'Location','southwest','FontSize',fs);
	
	title_name = sprintf('HARQ-%s : LDPC %d, Rate %.3f, max. iter %d, max. rounds %d qam mod %d mod app %d', combining_scheme, N, R, max_iter, max_rounds, qam_mod, mod_approx);
	title(title_name,'FontSize',fs);
	
	common_str = res_folder_all + sprintf('/BER_LDPC_%d_rate_%.3f_dec_iter_%d_max_rounds_%d_qm_%d_ma_%d', N,R, max_iter, max_rounds,qam_mod, mod_approx);
	
	fig_name = common_str + ".fig";
	savefig(fig_name);
	png_name = common_str + ".png";
	saveas(f,png_name);
	
	% set the figure properties
	figure('Renderer','painters','Position',[1000 400 800 500]);
	
	f = semilogy(SNRdB_vec,BLER_vec_pr_HARQ(end,:),'b-o');
	
	fs = 12;
	xlabel('SNR','FontSize',fs);
	ylabel('BLER','FontSize',fs);
	
	codeRate = R;
	leg_HARQ = sprintf('HARQ-%s BLER Rate %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
	
	legend(leg_HARQ, 'Location','southwest','FontSize',fs);
	
	title_name = sprintf('HARQ-%s : LDPC %d, Rate %.3f, max. iter %d, max. rounds %d qam mod %d mod app %d', combining_scheme, N, R, max_iter, max_rounds, qam_mod, mod_approx);
	title(title_name,'FontSize',fs);
	
	common_str = res_folder_all + sprintf('/BLER_LDPC_%d_rate_%.3f_dec_iter_%d_max_rounds_%d_qm_%d_ma_%d', N,R, max_iter, max_rounds,qam_mod, mod_approx);
	
	fig_name = common_str + ".fig";
	savefig(fig_name);
	png_name = common_str + ".png";
	saveas(f,png_name);
	
	%% Avg rounds

	% set the figure properties
	figure('Renderer','painters','Position',[1000 400 800 500]);
	
	f = semilogy(SNRdB_vec,Avg_rounds_HARQ,'b-o');
	
	fs = 12;
	xlabel('SNR','FontSize',fs);
	ylabel('Avg. rounds','FontSize',fs);
	
	leg_HARQ = sprintf('HARQ-%s AR FB Rate %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
	legend(leg_HARQ, 'Location','southwest','FontSize',fs);
	
	title_name = sprintf('HARQ-%s: LDPC %d, Rate %.3f, max. iter %d, max. rounds %d qam mod %d mod app %d', combining_scheme, N, R, max_iter, max_rounds, qam_mod, mod_approx);
	title(title_name,'FontSize',fs);
	
	common_str = res_folder_all + sprintf('/AR_LDPC_%d_rate_%.3f_dec_iter_%d_max_rounds_%d_qm_%d_ma_%d', N,R, max_iter, max_rounds,qam_mod, mod_approx);
	
	fig_name = common_str + ".fig";
	savefig(fig_name);
	png_name = common_str + ".png";
	saveas(f,png_name);
	
	% Store data    
	rate = R;
	data_file_name = res_folder_all + sprintf('/harq_data_LDPC_%d_rate_%.3f_dec_iter_%d_max_rounds_%d_qm_%d_ma_%d.mat', N,R, max_iter, max_rounds,qam_mod, mod_approx);
	save(data_file_name,'err_data','ar_data','snr_data','rate');
end

if (process_data_fb == 1)
    actualCodeRate_list(i_e) = R;
    %% Save and plot the results : all rounds
    err_data_all_thr(1,i_e,:,:) = BER_vec_pr_FB;
    err_data_all_thr(2,i_e,:,:) = BLER_vec_pr_FB;

    ar_data_all_thr(i_e,:) = Avg_rounds_FB;
    
    % Only current thr
    err_data = squeeze(err_data_all_thr(:,i_e,:,:));
    ar_data = squeeze(ar_data_all_thr(i_e,:));
    
    % BER
    % set the figure properties
    figure('Renderer','painters','Position',[1000 400 800 500]);
    
    f = semilogy(SNRdB_vec,BER_vec_pr_FB(end,:),'r-d');
    
    fs = 12;
    xlabel('SNR','FontSize',fs);
    ylabel('BER','FontSize',fs);
    
    codeRate = R;
    leg_FB = sprintf('FB-%s BER Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);
    
    legend(leg_FB, 'Location','southwest','FontSize',fs);
    
    title_name = sprintf('FB-%s : BER LDPC %d, mod. %s, Rate %.3f, max. iter %d, errthr %.3f, max. rounds %d qam mod %d mod app %d',combining_scheme, N, modulation, R, max_iter, err_thr, max_rounds, qam_mod, mod_approx);
    title(title_name,'FontSize',fs);
    
    common_str = res_folder_all + sprintf('/BER_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_max_rounds_%d_qm_%d_ma_%d', N,R, max_iter, err_thr, max_rounds,qam_mod, mod_approx);
    
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
    
    codeRate = R;
    leg_FB = sprintf('FB-%s BLER Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);
    
    legend(leg_FB, 'Location','southwest','FontSize',fs);
    
    title_name = sprintf('FB-%s : BLER LDPC %d, mod. %s, Rate %.3f, max. iter %d, errthr %.3f, max. rounds %d qam mod %d mod app %d',combining_scheme, N, modulation, R, max_iter, err_thr, max_rounds, qam_mod, mod_approx);
    title(title_name,'FontSize',fs);
    
    common_str = res_folder_all + sprintf('/BLER_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_max_rounds_%d_qm_%d_ma_%d', N,R, max_iter, err_thr, max_rounds,qam_mod, mod_approx);
    
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
    
    leg_FB = sprintf('FB-%s AR Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);
    legend(leg_FB, 'Location','southwest','FontSize',fs);
    
    title_name = sprintf('FB-%s : AR LDPC %d, mod. %s, Rate %.3f, max. iter %d, errthr %.3f, max. rounds %d qam mod %d mod app %d', combining_scheme, N, modulation, R, max_iter, err_thr, max_rounds, qam_mod, mod_approx);
    title(title_name,'FontSize',fs);
    
    common_str = res_folder_all + sprintf('/AR_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_max_rounds_%d_qm_%d_ma_%d', N,R, max_iter, err_thr, max_rounds,qam_mod, mod_approx);
    
    fig_name = common_str + ".fig";
    savefig(fig_name);
    png_name = common_str + ".png";
    saveas(f,png_name);
    
    % Store data 
    data_file_name = res_folder_all + sprintf('/fb_data_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_max_rounds_%d_qm_%d_ma_%d.mat', N,R, max_iter, err_thr, max_rounds,qam_mod, mod_approx);
    save(data_file_name,'err_data','ar_data','snr_data','err_thr');
end

%%
if (err_thr_ada_scheme == "opt")
    disp("err_thr_used : opt: ");
    disp(err_thr_ada_list);
else
    disp("err_thr_used : est: ");
    disp(err_thr_ada_list_est);
end
