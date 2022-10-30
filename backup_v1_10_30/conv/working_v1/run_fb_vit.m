% Fixing the range for reproducibility 
rng("default");

% Code parameters
targetCodeRate = 0.8;
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
mod_approx = 0;
comm_mod = 1;

err_thr_ada_scheme = "opt";
channel = "awgn";

% Error Stats
BER_vec_FB = zeros(1,num_SNRdB);
BLER_vec_FB = zeros(1,num_SNRdB);

BER_vec_pr_FB = zeros(max_rounds,num_SNRdB);
BLER_vec_pr_FB = zeros(max_rounds,num_SNRdB);

K = 400;

N = ceil(K/targetCodeRate);

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

nOut = 10;
nMiniFrames = 1000;

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
            [dataIn, rr_len] = conv_enc(data, targetCodeRate);
            
            % Symbol Modulation
            bpskModulator = comm.BPSKModulator;
            bpskModulator.PhaseOffset = pi/4;
            txSig = bpskModulator(dataIn);

            % Pass through channel
            rxSig = awgn(txSig,SNRdB);
            if (channel == "rayleigh")
                h = sqrt(rand(1)^2 + rand(1)^2);
                rxSig = awgn(h*txSig,SNRdB);
            end

            % Symbol demod
            if (mod_approx == 0)
                bpskDemodulator = comm.BPSKDemodulator; 
                bpskDemodulator.PhaseOffset = pi/4; 
                bpskDemodulator.DecisionMethod = 'Approximate log-likelihood ratio';
                rxLLR = -1*bpskDemodulator(rxSig);
            else
                bpskDemodulator = comm.BPSKDemodulator; 
                bpskDemodulator.PhaseOffset = pi/4; 
                bpskDemodulator.DecisionMethod = 'Hard decision';
                rxLLR = -1 + 2*bpskDemodulator(rxSig);
            end

            % Rate recovery and Decoding
            data_est = conv_dec(rxLLR, rr_len);

            % Check for error stats
            num_err = sum(data ~= double(data_est));
            num_err_FB = num_err;
            num_err_FB_per_round = zeros(max_rounds,1);
            num_err_FB_per_round(1) = num_err;
            
            % Start retransmission if 1st round failed
            
            if (num_err > 0 && max_rounds > 1)
                out = retransmit_func_FB_vit(channel,SNRdB,modulation,length(rxLLR),K,targetCodeRate,data,rxLLR,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,comm_mod,mod_approx,seed);
                num_ar_fb = num_ar_fb + out.Avg_rounds_FB;
                num_err_FB = out.num_err_FB;
                num_err_FB_per_round = out.num_err_vec;
            end

            BER_FB = BER_FB + num_err_FB;

            BLER_FB = BLER_FB + (num_err_FB > 0);

            BER_FB_per_round = BER_FB_per_round + num_err_FB_per_round;

            BLER_FB_per_round = BLER_FB_per_round + (num_err_FB_per_round > 0);
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
Avg_rounds_FB = Avg_rounds_FB/(i_on*nMiniFrames);

toc;

semilogy(SNRdB_vec,BLER_vec_FB); hold on;

%% 
disp("err_thr_used : ");
disp(err_thr_ada_list);