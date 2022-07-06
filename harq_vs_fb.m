% Fixing the range for reproducibility 
rng("default");

% Simulation settings : Defaults
if ~exist('global_settings','var')
    % Code parameters
    targetCodeRate = 0.8;
    nPRB = 8; % Vary this to change the code length
    max_iter = 6; % default is 8 in MATLAB
    err_thr = 0.05;
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
SNRdB_low = -8;
SNRdB_high = 10;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
num_snr = length(SNRdB_vec);

M = 2;
%     num_blocks = 10e1;

BER_vec = zeros(size(SNRdB_vec));
BLER_vec = zeros(size(SNRdB_vec));

BER_vec_ARQ = zeros(size(SNRdB_vec));
BLER_vec_ARQ = zeros(size(SNRdB_vec));

BER_vec_HARQ = zeros(size(SNRdB_vec));
BLER_vec_HARQ = zeros(size(SNRdB_vec));

BER_vec_FB = zeros(size(SNRdB_vec));
BLER_vec_FB = zeros(size(SNRdB_vec));

Avg_rounds_ARQ = zeros(size(SNRdB_vec));
Avg_rounds_HARQ = zeros(size(SNRdB_vec));
Avg_rounds_FB = zeros(size(SNRdB_vec));

% Retransmission params
max_rounds = 10;
fb_scheme = "HARQ";
harq_for_fb = 1;

% Compressor set up
p = 100*err_thr;
counts = [100-p p];

parpool_size = min(48,feature('numcores'));
if (parpool_size ~= 8)
    parpool(parpool_size);
end

tic;
parfor i_s = 1:length(SNRdB_vec)
    SNRdB = SNRdB_vec(i_s);
    noiseVar = 1./(10.^(SNRdB/10));

    BER = 0; BLER = 0;

    BER_ARQ = 0; BLER_ARQ = 0;

    BER_HARQ = 0; BLER_HARQ = 0;

    BER_FB = 0; BLER_FB = 0;

    % Print status
    fprintf('Status %0.2f %% done \n', round(i_s/num_snr*100));

    for i_n = 1:num_blocks
        
        % Update the counter 
        Avg_rounds_ARQ(i_s) = Avg_rounds_ARQ(i_s) + 1;
        Avg_rounds_HARQ(i_s) = Avg_rounds_HARQ(i_s) + 1; 
        Avg_rounds_FB(i_s) = Avg_rounds_FB(i_s) + 1;

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

        % Check for errors
        num_err = sum(mod(data+double(data_est),2));
        num_err_ARQ = num_err;
        num_err_HARQ = num_err;
        num_err_FB = num_err;
       
        if (num_err > 0)
            % Retransmission if failure
            % Always run ARQ for baseline
            for i_r = 1:max_rounds-1

                % Update the counter
                Avg_rounds_ARQ(i_s) = Avg_rounds_ARQ(i_s) + 1; 
                newRxSig = awgn(txSig,SNRdB);
                rxSig_ARQ = newRxSig
    
                % QAM Demod
                rxLLR_ARQ = qamdemod(rxSig_ARQ,M,'OutputType','LLR');
                        
                % Decoding and Rate recovery
                bgn = bgn_select(K,R);
                [data_est_ARQ, crc_chk_ARQ]  = nrldpc_dec(rxLLR_ARQ, R, modulation, K, max_iter, rv, nlayers, bgn);
        
                % Check for errors
                num_err_ARQ = sum(mod(data+double(data_est_ARQ),2));
                
                if (num_err_ARQ == 0)
                    break;
                end
            end

            % Always run HARQ for baseline
            prev_rxSig = rxSig;
            for i_r = 1:max_rounds-1

                % Update the counter
                Avg_rounds_HARQ(i_s) = Avg_rounds_HARQ(i_s) + 1; 
                newRxSig = awgn(txSig,SNRdB);
    
                % Chase combining
                prev_rxSig = [prev_rxSig newRxSig];
                rxSig_HARQ = mean(prev_rxSig,2);
    
                % QAM Demod
                rxLLR_HARQ = qamdemod(rxSig_HARQ,M,'OutputType','LLR');
                        
                % Decoding and Rate recovery
                bgn = bgn_select(K,R);
                [data_est_HARQ, crc_chk_HARQ]  = nrldpc_dec(rxLLR_HARQ, R, modulation, K, max_iter, rv, nlayers, bgn);
        
                % Check for errors
                num_err_HARQ = sum(mod(data+double(data_est_HARQ),2));
                
                if (num_err_HARQ == 0)
                    break;
                end
            end
            
            % Always run FB for comparison
            prev_rxSig = rxSig;
            prev_rxSig_FB = [];
            data_est_FB = data_est;
            decision_switch = 0;
            fb_count = 0; % number of times FB scheme was tried and failed
            fb_scheme = "HARQ";
            alg = 1;

            for i_r = 1:max_rounds-1
                num_err_FB = sum(mod(data+double(data_est_FB),2));

                % After every transmission, choose the retransmission scheme
                % Alg1 : HARQ over FB
                % Alg2 : FB over FB

                err_per = num_err_FB/K;
                % Once the error becomes sparse enough, stay on FB scheme
                if (decision_switch == 0)
                    if (err_per <= err_thr)
                        fb_scheme = "FB";
                        decision_switch = 1;
                    else
                        fb_scheme = "HARQ";
                    end
                end

                % Update the counter
                Avg_rounds_FB(i_s) = Avg_rounds_FB(i_s) + 1; 
                if (fb_scheme == "FB")
                    % Compress the error assuming free feedback
                    data_est_err = mod(data+double(data_est_FB),2);
                    err_seq = arithenco(data_est_err+1,counts);
    
                    assert(length(err_seq) < length(data_est_err),'Compression failed in FB round %d',i_r);
                    
                    targetErrCodeRate = length(err_seq)/N;

                    % Pick nPRB such that K_err >= length(err_seq)
                    tbs_err = nPRB_select(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate,length(err_seq));
                                        
                    K_err = tbs_err;
                    N_err = nPRB*NREPerPRB;
                    R_err = K_err/N_err;
    
                    bgn_err = bgn_select(K_err,R_err);
    
                    % Padding with zeroes
                    nz = K_err - length(err_seq);
                    assert(nz >= 0,'Zero padding failed in FB round %d',i_r);
                    err_seq_n = [err_seq; zeros(nz,1)];
                    
                    % Encoding and Rate matching
                    errDataInSeq = nrldpc_enc(err_seq_n, R_err, modulation, rv, bgn_err, nlayers);
    
                    % Reshape data into binary k-tuples, k = log2(M), and convert to
                    % integers
                    k = log2(M);
                    txData = reshape(errDataInSeq,length(errDataInSeq)/k,k);
                    txDataSym = bi2de(txData);
            
                    % QAM Modulation
                    txSig_FB = qammod(txDataSym,M);
                        
                    % Pass through AWGN channel
                    rxSig_FB = awgn(txSig_FB,SNRdB);

                    % Start storing the rxSig for HARQ
                    if harq_for_fb
                        prev_rxSig_FB = [prev_rxSig_FB rxSig_FB];
                        rxSig_FB = mean(prev_rxSig_FB,2);
                    end
            
                    % QAM Demod
                    rxLLR_FB = qamdemod(rxSig_FB,M,'OutputType','LLR');
                
                    % Decoding and Rate recovery
                    [err_seq_est, crc_chk_FB]  = nrldpc_dec(rxLLR_FB, R_err, modulation, K_err, max_iter, rv, nlayers, bgn_err);
    
                    % Decompress
                    err_deseq_est = arithdeco(double(err_seq_est),counts,length(data))-1;
                       
                    data_est_FB_new = mod(data_est_FB + err_deseq_est,2);
                    num_err_FB_new = sum(mod(data+double(data_est_FB_new),2));
                    % Correct the error only if CRC passes else round
                    % wasted
                    if (crc_chk_FB == 0 || num_err_FB_new == 0)
                        data_est_FB = mod(data_est_FB + err_deseq_est,2);
                        num_err_FB = sum(mod(data+double(data_est_FB),2));
                    end

                    % Empty the buffer if crc passes
                    if (crc_chk_FB == 0)
                        prev_rxSig_FB = []
                    end

                    if (num_err_FB_new == 0)
                        break;
                    end
                else
                    % Retransmit
                    newRxSig = awgn(txSig,SNRdB);
        
                    % Chase combining
                    prev_rxSig = [prev_rxSig newRxSig];
                    rxSig_FB = mean(prev_rxSig,2);
        
                    % QAM Demod
                    rxLLR_FB = qamdemod(rxSig_FB,M,'OutputType','LLR');
                            
                    % Decoding and Rate recovery
                    bgn = bgn_select(K,R);
                    [data_est_FB, crc_chk_FB] = nrldpc_dec(rxLLR_FB, R, modulation, K, max_iter, rv, nlayers, bgn);
            
                    % Check for errors
                    num_err_FB = sum(mod(data+double(data_est_FB),2));
    
                    if (num_err_FB == 0)
                        break;
                    end
                end

            end
        end

        BER = BER + num_err;
        BER_ARQ = BER_ARQ + num_err_HARQ;
        BER_HARQ = BER_HARQ + num_err_HARQ;
        BER_FB = BER_FB + num_err_FB;

        if (num_err > 0)
            BLER = BLER + 1;
        end
        if (num_err_ARQ > 0)
            BLER_ARQ = BLER_ARQ + 1;
        end

        if (num_err_HARQ > 0)
            BLER_HARQ = BLER_HARQ + 1;
        end

        if (num_err_FB > 0)
            BLER_FB = BLER_FB + 1;
        end
    end
    BER_vec(i_s) = BER/(K*num_blocks);
    BLER_vec(i_s) = BLER/(num_blocks);
    % ARQ
    BER_vec_ARQ(i_s) = BER_ARQ/(K*num_blocks);
    BLER_vec_ARQ(i_s) = BLER_ARQ/(num_blocks);
    % HARQ
    BER_vec_HARQ(i_s) = BER_HARQ/(K*num_blocks);
    BLER_vec_HARQ(i_s) = BLER_HARQ/(num_blocks);
    % FB
    BER_vec_FB(i_s) = BER_FB/(K*num_blocks);
    BLER_vec_FB(i_s) = BLER_FB/(num_blocks);

    % Avg rounds per transmission
    Avg_rounds_ARQ(i_s) = Avg_rounds_ARQ(i_s)/num_blocks;
    Avg_rounds_HARQ(i_s) = Avg_rounds_HARQ(i_s)/num_blocks;
    Avg_rounds_FB(i_s) = Avg_rounds_FB(i_s)/num_blocks;
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
semilogy(SNRdB_vec,BLER_vec_FB);
xlabel('SNR');
ylabel('BLER');
leg_ARQ = sprintf('ARQ max. %d rounds',max_rounds);
leg_HARQ = sprintf('HARQ max. %d rounds',max_rounds);
leg_FB = sprintf('Feedback max. %d rounds',max_rounds);
legend('No HARQ', leg_ARQ, leg_HARQ, leg_FB);
title_name = sprintf('BLER HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter err thr %.2f',N,K,R, max_iter, err_thr);
title(title_name);
filename_BLER_fig = sprintf('results/BLER_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f.fig',N,R, max_iter, err_thr);
filename_BLER_png = sprintf('results/BLER_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f.png',N,R, max_iter, err_thr);
savefig(filename_BLER_fig);
saveas(f,filename_BLER_png);

figure('Renderer','painters','Position',[1000 400 800 500]);
f = semilogy(SNRdB_vec,Avg_rounds_ARQ);
hold on;
semilogy(SNRdB_vec,Avg_rounds_HARQ);
semilogy(SNRdB_vec,Avg_rounds_FB);
xlabel('SNR');
ylabel('Avg rounds');
leg_ARQ = sprintf('ARQ max. %d rounds',max_rounds);
leg_HARQ = sprintf('HARQ max. %d rounds',max_rounds);
leg_FB = sprintf('Feedback max. %d rounds',max_rounds);
legend(leg_ARQ, leg_HARQ, leg_FB);
title_name = sprintf('HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter %.2f',N,K,R, max_iter, err_thr);
title(title_name);
filename_AR_fig = sprintf('results/AR_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f.fig',N,R, max_iter, err_thr);
filename_AR_png = sprintf('results/AR_LDPC_%d_rate_%.2f_dec_iter_%d_err_thr_%.2f.png',N,R, max_iter, err_thr);
savefig(filename_AR_fig);
saveas(f,filename_AR_png);

%%
function [output] = nrldpc_enc(data_in, rate, modulation, rv, bgn, nlayers)
    in = int8(data_in);
    tbIn = nrCRCEncode(in,  '16');
    cbsIn = nrCodeBlockSegmentLDPC(tbIn, bgn);
    enc = nrLDPCEncode(cbsIn, bgn);
    outlen = ceil(length(data_in)/rate);
    output = nrRateMatchLDPC(enc,outlen,rv,modulation,nlayers);
    output = double(output);
end

function [output, crc] = nrldpc_dec(chOut, rate, modulation, block_len, max_iter, rv, nlayers, bgn)
    CRC_L = 16;

    raterec = nrRateRecoverLDPC(chOut, block_len, rate, rv, modulation, nlayers);

    decBits = nrLDPCDecode(raterec, bgn, max_iter);

    % Code block desegmentation and CRC decoding
    [blk, ~] = nrCodeBlockDesegmentLDPC(decBits, bgn, block_len+CRC_L);

    % Transport block CRC decoding
    [output, crc] = nrCRCDecode(blk,  '16');
    output = double(output);
end
