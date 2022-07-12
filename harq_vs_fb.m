% Fixing the range for reproducibility 
rng("default");

% Simulation settings : Defaults
if ~exist('global_settings','var')
    % Code parameters
    targetCodeRate = 0.8;
    nPRB = 20; % Vary this to change the code length
    max_iter = 6; % default is 8 in MATLAB
    min_ber = 1e-5;
    SNRdB_low = -8;%10;
    nFrames = 10e1;
    err_thr_ada = 0;
    err_thr = 0.04;
end

% use nrTBS to get K,N
modulation = 'pi/2-BPSK';
nlayers = 1;

NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice
tbs = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetCodeRate);

% Use generated code params
K = tbs;
N = nPRB*NREPerPRB;
R = K/N;

% base graph selection based on rate and K
bgn = bgn_select(K,R);

rv = 0;

% Simulation params
SNRdB_step = 0.2;
SNRdB_high = SNRdB_low+18;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
num_snr = length(SNRdB_vec);

M = 2;

BER_vec = zeros(size(SNRdB_vec));
BLER_vec = zeros(size(SNRdB_vec));

BER_vec_ARQ = zeros(size(SNRdB_vec));
BLER_vec_ARQ = zeros(size(SNRdB_vec));

BER_vec_HARQ = zeros(size(SNRdB_vec));
BLER_vec_HARQ = zeros(size(SNRdB_vec));

BER_vec_FB_HARQ = zeros(size(SNRdB_vec));
BLER_vec_FB_HARQ = zeros(size(SNRdB_vec));

BER_vec_FB_FB = zeros(size(SNRdB_vec));
BLER_vec_FB_FB = zeros(size(SNRdB_vec));

Avg_rounds_ARQ = zeros(size(SNRdB_vec));
Avg_rounds_HARQ = zeros(size(SNRdB_vec));
Avg_rounds_FB_HARQ = zeros(size(SNRdB_vec));
Avg_rounds_FB_FB = zeros(size(SNRdB_vec));

% Retransmission params
max_rounds = 10;
fb_scheme = "HARQ";
harq_for_fb = 1;

parpool_size = min(64,feature('numcores'));
if (parpool_size ~= 8)
    parpool(parpool_size);
end

if (err_thr_ada)
    for id = 1:length(SNRdB_vec)
        err_thr_ada_list(id) = err_thr_select(N,R,round(SNRdB_vec(id),4),min_ber);
    end
end

err_thr_ada_list = err_thr*ones(size(SNRdB_vec));
tic;
parfor i_s = 1:length(SNRdB_vec)
    SNRdB = round(SNRdB_vec(i_s),4);

    noiseVar = 1./(10.^(SNRdB/10));

    err_thr = err_thr_ada_list(i_s);

    % Compressor set up
    p = round(100*err_thr);
    counts = [100-p p];

    BER = 0; BLER = 0;

    BER_ARQ = 0; BLER_ARQ = 0;

    BER_HARQ = 0; BLER_HARQ = 0;

    BER_FB_HARQ = 0; BLER_FB_HARQ = 0;

    BER_FB_FB = 0; BLER_FB_FB = 0;

    % Print status
    fprintf('Status %0.2f %% done \n', round(i_s/num_snr*100));

    for i_n = 1:nFrames
        
        % Update the counter 
        Avg_rounds_ARQ(i_s) = Avg_rounds_ARQ(i_s) + 1;
        Avg_rounds_HARQ(i_s) = Avg_rounds_HARQ(i_s) + 1; 
        Avg_rounds_FB_HARQ(i_s) = Avg_rounds_FB_HARQ(i_s) + 1;
        Avg_rounds_FB_FB(i_s) = Avg_rounds_FB_FB(i_s) + 1;

        % Generate random message
        data = randi([0 1], K, 1);
    
        % Encoding and Rate matching
        bgn = bgn_select(K,R);
        dataIn = nrldpc_enc(data, R, modulation, rv, bgn, nlayers);

        % Reshape data into binary k-tuples, k = log2(M), and convert to
        % integers
        k = log2(M);
        txData = reshape(dataIn,length(dataIn)/k,k);
        txDataSym = bi2de(txData);

        % QAM Modulation
        txSig = qammod(txDataSym,M);
            
        % Pass through AWGN channel
        rxSig = awgn(txSig,SNRdB);

        % QAM Demod
        rxLLR = qamdemod(rxSig,M,'OutputType','LLR');
                
        % Decoding and Rate recovery
        bgn = bgn_select(K,R);
        [data_est, crc_chk] = nrldpc_dec(rxLLR, R, modulation, K, max_iter, rv, nlayers, bgn);

        % Check for error stats
        num_err = sum(mod(data+double(data_est),2));
        num_err_ARQ = num_err;
        num_err_HARQ = num_err;
        num_err_FB_HARQ = num_err;
        num_err_FB_FB = num_err;
        
        % Start retransmission if 1st round failed
        if (crc_chk > 0)

            out = retransmit_func_ARQ(SNRdB,modulation,max_iter,rv,nlayers,M,K,R,data,txSig,max_rounds);
            Avg_rounds_ARQ(i_s) = Avg_rounds_ARQ(i_s) + out.Avg_rounds_ARQ;
            num_err_ARQ = out.num_err_ARQ;

            out = retransmit_func_HARQ(SNRdB,modulation,max_iter,rv,nlayers,M,K,R,data,txSig,rxSig,max_rounds);
            Avg_rounds_HARQ(i_s) = Avg_rounds_HARQ(i_s) + out.Avg_rounds_HARQ;
            num_err_HARQ = out.num_err_HARQ;
            
            out = retransmit_func_FB_HARQ(SNRdB,modulation,max_iter,rv,nlayers,nPRB,NREPerPRB,M,N,K,R,data,txSig,rxSig,data_est,err_thr,max_rounds,counts,num_err);
            Avg_rounds_FB_HARQ(i_s) = Avg_rounds_FB_HARQ(i_s) + out.Avg_rounds_FB_HARQ;
            num_err_FB_HARQ = out.num_err_FB_HARQ; 

            out = retransmit_func_FB_FB(SNRdB,modulation,max_iter,rv,nlayers,nPRB,NREPerPRB,M,N,K,R,data,txSig,rxSig,data_est,err_thr,max_rounds,counts,num_err);
            Avg_rounds_FB_FB(i_s) = Avg_rounds_FB_FB(i_s) + out.Avg_rounds_FB_FB;
            num_err_FB_FB = out.num_err_FB_FB;
        end

        BER = BER + num_err;
        BER_ARQ = BER_ARQ + num_err_HARQ;
        BER_HARQ = BER_HARQ + num_err_HARQ;
        BER_FB_HARQ = BER_FB_HARQ + num_err_FB_HARQ;
        BER_FB_FB = BER_FB_FB + num_err_FB_FB;

        if (num_err > 0)
            BLER = BLER + 1;
        end
        if (num_err_ARQ > 0)
            BLER_ARQ = BLER_ARQ + 1;
        end

        if (num_err_HARQ > 0)
            BLER_HARQ = BLER_HARQ + 1;
        end

        if (num_err_FB_HARQ > 0)
            BLER_FB_HARQ = BLER_FB_HARQ + 1;
        end

        if (num_err_FB_FB > 0)
            BLER_FB_FB = BLER_FB_FB + 1;
        end
    end
    BER_vec(i_s) = BER/(K*nFrames);
    BLER_vec(i_s) = BLER/(nFrames);
    % ARQ
    BER_vec_ARQ(i_s) = BER_ARQ/(K*nFrames);
    BLER_vec_ARQ(i_s) = BLER_ARQ/(nFrames);
    % HARQ
    BER_vec_HARQ(i_s) = BER_HARQ/(K*nFrames);
    BLER_vec_HARQ(i_s) = BLER_HARQ/(nFrames);
    % FB HARQ
    BER_vec_FB_HARQ(i_s) = BER_FB_HARQ/(K*nFrames);
    BLER_vec_FB_HARQ(i_s) = BLER_FB_HARQ/(nFrames);
    % FB FB
    BER_vec_FB_FB(i_s) = BER_FB_FB/(K*nFrames);
    BLER_vec_FB_FB(i_s) = BLER_FB_FB/(nFrames);

    % Avg rounds per transmission
    Avg_rounds_ARQ(i_s) = Avg_rounds_ARQ(i_s)/nFrames;
    Avg_rounds_HARQ(i_s) = Avg_rounds_HARQ(i_s)/nFrames;
    Avg_rounds_FB_HARQ(i_s) = Avg_rounds_FB_HARQ(i_s)/nFrames;
    Avg_rounds_FB_FB(i_s) = Avg_rounds_FB_FB(i_s)/nFrames;
end
toc;
if (parpool_size ~= 8)
    delete(gcp('nocreate'));
end

%%
% Plots
figure('Renderer','painters','Position',[100 400 800 500]);
f = semilogy(SNRdB_vec,BLER_vec);
hold on;
semilogy(SNRdB_vec,BLER_vec_ARQ);
semilogy(SNRdB_vec,BLER_vec_HARQ);
semilogy(SNRdB_vec,BLER_vec_FB_HARQ);
semilogy(SNRdB_vec,BLER_vec_FB_FB);
xlabel('SNR');
ylabel('BLER');
leg_ARQ = sprintf('ARQ max. %d rounds',max_rounds);
leg_HARQ = sprintf('HARQ max. %d rounds',max_rounds);
leg_FB_HARQ = sprintf('FB+HARQ max. %d rounds',max_rounds);
leg_FB_FB = sprintf('FB+FB max. %d rounds',max_rounds);
legend('No HARQ', leg_ARQ, leg_HARQ, leg_FB_HARQ, leg_FB_FB);

if (err_thr_ada)
    title_name = sprintf('BLER HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter, adaptive err thr %.2f',N,K,R, max_iter, min_ber);
    title(title_name);
    filename_BLER_fig = sprintf('results/BLER_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_ada_%.2f_max_round_%d.fig',N,R, max_iter, min_ber, max_rounds);
    filename_BLER_png = sprintf('results/BLER_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_ada_%.2f_max_round_%d.png',N,R, max_iter, min_ber, max_rounds);
else
    title_name = sprintf('BLER HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter err thr %.2f',N,K,R, max_iter, err_thr);
    title(title_name);
    filename_BLER_fig = sprintf('results/BLER_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f_max_round_%d.fig',N,R, max_iter, err_thr, max_rounds);
    filename_BLER_png = sprintf('results/BLER_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f_max_round_%d.png',N,R, max_iter, err_thr, max_rounds);
end
savefig(filename_BLER_fig);
saveas(f,filename_BLER_png);

figure('Renderer','painters','Position',[1000 400 800 500]);
f = plot(SNRdB_vec,ones(size(SNRdB_vec)));
hold on;
plot(SNRdB_vec,Avg_rounds_ARQ);
plot(SNRdB_vec,Avg_rounds_HARQ);
plot(SNRdB_vec,Avg_rounds_FB_HARQ);
plot(SNRdB_vec,Avg_rounds_FB_FB);
xlabel('SNR');
ylabel('Avg rounds');
leg_ARQ = sprintf('ARQ max. %d rounds',max_rounds);
leg_HARQ = sprintf('HARQ max. %d rounds',max_rounds);
leg_FB_HARQ = sprintf('FB+HARQ max. %d rounds',max_rounds);
leg_FB_FB = sprintf('FB+FB max. %d rounds',max_rounds);
legend('No HARQ',leg_ARQ, leg_HARQ, leg_FB_HARQ, leg_FB_FB);

if (err_thr_ada)
    title_name = sprintf('Avg. rounds HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter, adaptive err thr %0.2f',N,K,R, max_iter, min_ber);
    title(title_name);
    filename_AR_fig = sprintf('results/AR_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_ada_%0.2f_max_round_%d.fig',N,R, max_iter, min_ber, max_rounds);
    filename_AR_png = sprintf('results/AR_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_ada_%0.2f_max_round_%d.png',N,R, max_iter, min_ber, max_rounds);
else
    title_name = sprintf('Avg. rounds vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter err thr %.2f',N,K,R, max_iter, err_thr);
    title(title_name);
    filename_AR_fig = sprintf('results/AR_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f_max_round_%d.fig',N,R, max_iter, err_thr, max_rounds);
    filename_AR_png = sprintf('results/AR_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f_max_round_%d.png',N,R, max_iter, err_thr, max_rounds);
end
savefig(filename_AR_fig);
saveas(f,filename_AR_png);


%% 
disp("err_thr_used : ");
disp(err_thr_ada_list);