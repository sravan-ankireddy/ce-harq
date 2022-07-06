% Fixing the range for reproducibility 
rng("default");

% Code parameters
targetCodeRate = 0.9;

% use nrTBS to get K,N
modulation = 'pi/2-BPSK';
nlayers = 1;
nPRB = 8; % Vary this to change the code rate
NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice
tbs = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetCodeRate);

% Use generated code params
K = tbs;
N = nPRB*NREPerPRB;
R = K/N;

% base graph selection based on rate and K
bgn = bgn_select(K,R);

rv = 0;
max_iter = 6; % default is 8 in MATLAB
err_thr = 0.02;

% Simulation params
SNRdB_step = 1;%0.25;
SNRdB_low = -2;%-5;
SNRdB_high = 5;%5;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;

M = 2;
num_blocks = 10e3;

BER_vec = zeros(size(SNRdB_vec));
BLER_vec = zeros(size(SNRdB_vec));

BER_vec_HARQ = zeros(size(SNRdB_vec));
BLER_vec_HARQ = zeros(size(SNRdB_vec));

BER_vec_FB = zeros(size(SNRdB_vec));
BLER_vec_FB = zeros(size(SNRdB_vec));

Avg_rounds_HARQ = zeros(size(SNRdB_vec));
Avg_rounds_FB = zeros(size(SNRdB_vec));

% Retransmission params
max_rounds = 10;
fb_scheme = "HARQ";

% Compressor set up
p = 100*err_thr;
counts = [100-p p];
tic;
parfor i_s = 1:length(SNRdB_vec)
    SNRdB = SNRdB_vec(i_s);
    noiseVar = 1./(10.^(SNRdB/10));
    BER = 0;
    BLER = 0;
    BER_HARQ = 0;e
    BLER_HARQ = 0;
    BER_FB = 0;
    BLER_FB = 0;

    % Print status
    fprintf('Status %0.2f %% done \n', round(i_s/length(SNRdB_vec)*100));

    for i_n = 1:num_blocks
        
        % Update the counter 
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
        num_err_HARQ = num_err;
        num_err_FB = num_err;
       
        if (num_err > 0)
            % Retransmission if failure
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
            data_est_FB = data_est;
            decision_switch = 0;
            fb_count = 0; % number of times FB scheme was tried and failed
            max_error_tx = 10;
            fb_scheme = "HARQ";
            alg = 2;
            chaseError = true;
            num_err_last_FB = intmax; % so that error compression is done in 1st feedback
            useNew = true; % signal to decoder whether to use last round or not - if CRC failed
            for i_r = 1:max_rounds-1
                num_err_FB = sum(mod(data+double(data_est_FB),2));

                % After every transmission, choose the retransmission scheme
                % FIX ME : 
                % Alg1 : Once switched from HARQ to FB, no swicthing back
                % Alg2 : After trying fb for 2 rounds, go back to HARQ is
                % it fails

                err_per = num_err_FB/K;
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

                        % If previous transmission had lesser number of
                        % errors, reuse it instead of compressing new
                        % error vector
                    if (chaseError && num_err_FB >= num_err_last_FB)
                        % skip all encoding, etc. Send symbol computed
                        % in previous round
                        new_rxSig_FB = awgn(txSig,SNRdB);
                        prev_rxSig_FB = [prev_rxSig_FB new_rxSig_FB];
                        rxSig_FB = mean(prev_rxSig_FB, 2);
                        useNew = false; % too many errors in previous transmission
                    else
                        err_seq = arithenco(data_est_err+1,counts);
        
                        assert(length(err_seq) < length(data_est_err),'Compression failed in FB round %d',i_r);
                        
                        targetErrCodeRate = length(err_seq)/N;
    
                        % Pick nPRB such that K_err >= length(err_seq)
                        tbs_err = nPRB_select(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate,length(err_seq));
                                            
                        K_err = tbs_err;
                        N_err = nPRB*NREPerPRB;
                        R_err = K_err/N_err;
        
                        bgn_err = bgn_select(K_err,N_err);
        
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
                        txSig = qammod(txDataSym,M);
                            
                        % Pass through AWGN channel
                        rxSig_FB = awgn(txSig,SNRdB);

                        if chaseError
                            prev_rxSig_FB = rxSig_FB;
                            useNew = true;
                        end
                    end
            
                    
                    % QAM Demod
                    rxLLR_FB = qamdemod(rxSig_FB,M,'OutputType','LLR');
                
                    % Decoding and Rate recovery
                    [err_seq_est, crc_chk_FB]  = nrldpc_dec(rxLLR_FB, R_err, modulation, K_err, max_iter, rv, nlayers, bgn_err);
    
                    % Decompress
                    err_deseq_est = arithdeco(double(err_seq_est),counts,length(data))-1;
                       
                    data_est_FB_new = mod(data_est_FB + err_deseq_est,2);
                    num_err_FB_new = sum(mod(data+double(data_est_FB),2));


                    if chaseError
                        % Correct even if CRC fails. If too many errors,
                        % can fall back to previous estimate.

                        if useNew % previous correction resulted in lesser errors
                            num_err_last_FB = num_err_FB;
                            data_est_FB = mod(data_est_FB + err_deseq_est,2);
                            num_err_FB = sum(mod(data+double(data_est_FB),2));
                            old_data_est_FB = data_est_FB;
                            
                        else % previous correction resulted in more errors: throw estimate of previous round
                            % => num_err_last_FB doesn't change
                            data_est_FB = mod(old_data_est_FB + err_deseq_est,2);
                            num_err_FB = sum(mod(data+double(old_data_est_FB),2));
                        end
                    else
                        % Correct the error only if CRC passes else round
                        % wasted
                        if (crc_chk_FB == 0 || num_err_FB_new == 0)
                            data_est_FB = mod(data_est_FB + err_deseq_est,2);
                            num_err_FB = sum(mod(data+double(data_est_FB),2));
                        end
                    end
                    
                    if (crc_chk_FB)
                        fb_count = fb_count + 1;
                        if (fb_count == max_error_tx)
                            fb_scheme = "HARQ";
                        end
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
        BER_HARQ = BER_HARQ + num_err_HARQ;
        BER_FB = BER_FB + num_err_FB;

        if (num_err > 0)
            BLER = BLER + 1;
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
    % HARQ
    BER_vec_HARQ(i_s) = BER_HARQ/(K*num_blocks);
    BLER_vec_HARQ(i_s) = BLER_HARQ/(num_blocks);
    % FB
    BER_vec_FB(i_s) = BER_FB/(K*num_blocks);
    BLER_vec_FB(i_s) = BLER_FB/(num_blocks);

    % Avg rounds per transmission
    Avg_rounds_HARQ(i_s) = Avg_rounds_HARQ(i_s)/num_blocks;
    Avg_rounds_FB(i_s) = Avg_rounds_FB(i_s)/num_blocks;
end
toc;

%%
figure(1)
% subplot(2,1,1);
% semilogy(SNRdB_vec,BER_vec);
% hold on;
% semilogy(SNRdB_vec,BER_vec_HARQ);
% semilogy(SNRdB_vec,BER_vec_FB);
% leg_HARQ = sprintf('HARQ max. %d rounds',max_rounds);
% leg_FB = sprintf('Feedback max. %d rounds',max_rounds);
% legend('No HARQ', leg_HARQ, leg_FB);
% title_name = sprintf('BER HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter',N,K,R, max_iter);
% title(title_name);

% subplot(2,1,2)
semilogy(SNRdB_vec,BLER_vec);
hold on;
semilogy(SNRdB_vec,BLER_vec_HARQ);
semilogy(SNRdB_vec,BLER_vec_FB);
xlabel('SNR');
ylabel('BLER');
leg_HARQ = sprintf('HARQ max. %d rounds',max_rounds);
leg_FB = sprintf('Feedback max. %d rounds',max_rounds);
legend('No HARQ', leg_HARQ, leg_FB);
grid on
title_name = sprintf('BLER HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter',N,K,R, max_iter);
title(title_name);

figure(2)
semilogy(SNRdB_vec,Avg_rounds_HARQ);
hold on;
semilogy(SNRdB_vec,Avg_rounds_FB);
xlabel('SNR');
ylabel('Avg rounds');
leg_HARQ = sprintf('HARQ max. %d rounds',max_rounds);
leg_FB = sprintf('Feedback max. %d rounds',max_rounds);
legend(leg_HARQ, leg_FB);
grid on
title_name = sprintf('HARQ vs Feedback : LDPC (%d,%d), Rate %.2f, %d decoding iter',N,K,R, max_iter);
title(title_name);
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
