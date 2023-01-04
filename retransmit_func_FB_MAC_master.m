function out = retransmit_func_FB_MAC_master(channel,SNRdB,modulation,N,K,R,code_type,feedback_mode,combining_scheme,rvSeq,ncb,Nref,max_iter,nlayers,dec_type,data,rxLLR,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,comm_mod,mod_approx,seed)
    
    rng(seed);
    int_state = seed;
    base_rate = R;
    M = bits_per_symbol(modulation);
    Avg_rounds_FB = 0;
    rxLLR_HARQ_buffer = rxLLR;
    rxLLR_FB_buffer = [];
    rxLLR_FB_buffer_outer = [];
    data_est_FB = data_est;
    data_est_FB_prev = data_est;
    decision_switch = 0;
    fb_scheme = "HARQ";
    num_err_FB_prev = num_err;
    num_err_vec = zeros(max_rounds,1);
    num_err_vec(1) = num_err;

    noiseVar = 1./(10.^(SNRdB/10));
    rep_coding = 0;
    tc_rp = 0.2;

    comp_rates = 1./(2:12);
    minR = comp_rates(end);
    rc = 0;

    % FIX ME
    rv = 0;
    rvSeq_ind = 1;

    min_bler = 1e-4;
    harq_data_path = sprintf('bler_data/awgn/%d/%s/harq/%s/100000/harq_data_Conv_%d_rate_0.833_rate_0.083_max_rounds_%d.mat', N, dec_type, modulation, N, max_rounds);
    acomp_table_path = sprintf('lut_data/acomp_%d_ns_100000.mat', K);
    if (err_thr_ada_scheme == "est")
        harq_data = load(harq_data_path);
        acomp_table = load(acom_table_path);
    end
    prev_reset_round = 0;%max_rounds;
    for i_r = 1:max_rounds-1

        num_err_FB = sum(data ~= double(data_est_FB));
        err_per = num_err_FB/K;

        % Adaptive err thr
        if (err_thr_ada_scheme == "est")
            % err_thr = err_thr_ada_list_est((max_rounds - i_r),i_s);
            if (fb_scheme == "FB")
                err_thr = err_thr_select(harq_data,acomp_table,targetErrCodeRate,SNRdB,max_rounds - prev_reset_round,max_rounds - i_r,min_bler);
            else
                err_thr = err_thr_select(harq_data,acomp_table,R,SNRdB,max_rounds,max_rounds - i_r,min_bler);
            end
            % err_thr_ada_list_est((max_rounds - i_r),i_s) = err_thr;
        end

        % Once the error becomes sparse enough, stay on FB scheme
        % Within the FB scheme, you can choose FB/HARQ based on error
        % vector
        if (decision_switch == 0)
            if (err_per <= err_thr && err_thr > 0)
                % disp("found sparse err ")
                % Also check if the error is compressible
                data_est_err_temp = mod(data+double(data_est_FB),2);
                % iterative counts
                p_cur = max(1,round(100 * sum(data_est_err_temp)/length(data_est_err_temp)));
                counts_temp = [100-p_cur p_cur];
                err_seq_temp = arithenco(data_est_err_temp+1,counts_temp);
                % err_seq_temp = arithenco(data_est_err_temp+1,counts);
                if (length(err_seq_temp) < length(data_est_err_temp))
                    fb_scheme = "FB";
                    decision_switch = 1;
                    counts = counts_temp;
                    % disp(" switched to feedback ")
                else
                    fb_scheme = "HARQ";
                end
            else
                fb_scheme = "HARQ";
            end
        end

        % Update the counter
        Avg_rounds_FB = Avg_rounds_FB + 1; 
        
        if (fb_scheme == "FB")
            % Compress the error assuming free feedback
            data_est_err = mod(data+double(data_est_FB),2);

            % Choose HARQ/FB based on error vector
            num_err_FB_cur = sum(data_est_err);
            if (num_err_FB_cur > num_err_FB_prev)
                % No. of errors increased from prev. round; use HARQ on FB
                inner_scheme = "HARQ";

                % Reset latest est to prev round est
                data_est_FB = data_est_FB_prev;
                num_err_FB = num_err_FB_prev;
            else
                inner_scheme = "FB";
                % Empty the chase buffer if prev round reduced the errors

                rxLLR_FB_buffer = [];
                rxLLR_FB_buffer_outer = [];
                % Reset rvSeq start
                rvSeq_ind = 1;
            end

            if (inner_scheme == "FB")
                % disp("inner feedback scheme ")
                
                % iterative counts
                p_cur = max(1,round(100 * sum(data_est_err)/length(data_est_err)));
                counts = [100-p_cur p_cur];
                err_seq = arithenco(data_est_err+1,counts);
    
                assert(length(err_seq) < length(data_est_err),'Compression failed in FB round %d',i_r);

                % find the available rate for inner code
                if (code_type == "Conv")
                    % PHY-MAC mode
                    if (feedback_mode == "MAC")
                        len_available = K;
                    else
                        len_available = N;
                    end
                    % find the smallest rate for compression : coding will be only for comp -> K, not N
                    targetErrCodeRate = (length(err_seq)+3)/len_available;
                    if (targetErrCodeRate >= minR)
                        % get the smallest rate possible
                        [~, min_ind] = min(abs(targetErrCodeRate - comp_rates));
                        if (comp_rates(min_ind) > targetErrCodeRate)
                            comp_rate = comp_rates(min_ind);
                        else
                            while (comp_rates(min_ind) < targetErrCodeRate)
                                min_ind = min_ind - 1;
                            end
                            comp_rate = comp_rates(min_ind);
                        end
                    else
                        comp_rate = comp_rates(end);
                        min_ind = length(comp_rates) - 1;
                    end
                    if (feedback_mode == "MAC")
                        % encode to length <= K
                        comp_code_inner = conv_enc(err_seq, comp_rate);

                        % interleaver between inner and outer code
                        comp_code_inner = randintrlv(comp_code_inner,int_state);

                        % append zeros to match length to K
                        nz = K - length(comp_code_inner);
                        comp_code_inner = [comp_code_inner; zeros(nz,1)];

                        % perform outer coding with rate K/N
                        comp_code_outer = conv_enc(comp_code_inner, base_rate);

                        errDataInSeq = comp_code_outer;
                    else
                        % Directly encode to length <= N
                        errDataInSeq = conv_enc(err_seq, comp_rate);
                    end
                else
                    % FIX ME : only PHY now for LDPC
                    % find the smallest rate for compression : coding will be only for comp -> K, not N
                    targetErrCodeRate = length(err_seq)/N;
                    k = bits_per_symbol(modulation);
                    M = 2^k;
                    % PRB settings
                    tarCodeLen = 960;
                    nlayers = 1;
                    NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice
                    nPRB = round(tarCodeLen/(k*NREPerPRB)); % Vary this to change the code length

                    % Pick nPRB such that K_err >= length(err_seq)
                    if (modulation == 'BPSK')
                        nPRB_out = nPRB_select('pi/2-BPSK',nlayers,nPRB,NREPerPRB,targetErrCodeRate,length(err_seq));
                    else
                        nPRB_out = nPRB_select(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate,length(err_seq));
                    end
                    tbs_err = nPRB_out.tbs_err;
                    nPRB = nPRB_out.nPRB;
                    
                    K_err = tbs_err; %no. bits in transportBlock
                    N_err = nPRB*NREPerPRB*k;
                    R_err = K_err/N_err;
        
                    bgn_err = bgn_select(K_err,R_err);
        
                    % Padding with zeroes
                    nz = K_err - length(err_seq);
                    err_seq_n = [err_seq; zeros(nz,1)];
                    
                    % Encoding and Rate matching
                    errDataInSeq = nrldpc_enc(err_seq_n, R_err, modulation, rv, bgn_err, nlayers);
                end

                % Symbol Modulation
                if (modulation == "BPSK")
                    bpskModulator = comm.BPSKModulator;
                    bpskModulator.PhaseOffset = pi/4;
                    txSig_FB = bpskModulator(errDataInSeq);
                else
                    txSig_FB = qammod(errDataInSeq,M,'InputType','bit','UnitAveragePower',true);
                end
                    
                % Pass through channel
                if (channel == "awgn")
                    rxSig_FB = awgn(txSig_FB,SNRdB,'measured');
                elseif (channel == "rayleigh")
                    h = sqrt(rand(1)^2 + rand(1)^2);
                    rxSig_FB = awgn(h*txSig_FB,SNRdB,'measured');
                end
            else
                %% FIX ME : Sticking to CC for now
 
                % Symbol Modulation
                if (modulation == "BPSK")
                    bpskModulator = comm.BPSKModulator;
                    bpskModulator.PhaseOffset = pi/4;
                    txSig_FB = bpskModulator(errDataInSeq);
                else
                    txSig_FB = qammod(errDataInSeq,M,'InputType','bit','UnitAveragePower',true);
                end

                % Pass through channel
                if (channel == "awgn")
                    rxSig_FB = awgn(txSig_FB,SNRdB,'measured');
                elseif (channel == "rayleigh")
                    h = sqrt(rand(1)^2 + rand(1)^2);
                    rxSig_FB = awgn(h*txSig_FB,SNRdB,'measured');
                end
            end
    
            % Symbol demod
            if (modulation == "BPSK")
                if (dec_type == "hard")
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Hard decision';
                    rxLLR_FB = bpskDemodulator(rxSig_FB);
                else
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Approximate log-likelihood ratio';
                    rxLLR_FB = bpskDemodulator(rxSig_FB);
                end
            else
                if (dec_type == "hard")
                    rxLLR_FB = qamdemod(rxSig_FB, M, OutputType='bit', UnitAveragePower=true);
                else
                    rxLLR_FB = qamdemod(rxSig_FB, M, OutputType='approxllr', UnitAveragePower=true, NoiseVariance=noiseVar);
                end
            end

            % Chase Combining
            if (dec_type == "hard")
                rxLLR_FB_buffer = [rxLLR_FB_buffer rxLLR_FB];
                rxLLR_FB = round(mean(rxLLR_FB_buffer, 2));
            else
                rxLLR_FB_buffer = [rxLLR_FB_buffer rxLLR_FB];
                rxLLR_FB = sum(rxLLR_FB_buffer, 2);
            end

            % Decoding
            if (code_type == "Conv")
                if (feedback_mode == "MAC")
                    % first decode outer code
                    outer_err_seq_est = conv_dec(rxLLR_FB,base_rate,dec_type);

                    % remove the zero padding
                    outer_err_seq_est = outer_err_seq_est(1:end-nz);
                    
                    % deinterleaver between inner and outer code
                    outer_err_seq_est = randdeintrlv(outer_err_seq_est,int_state);
                    
                    % next decode inner : map to 1 - 2*c for soft
                    if (dec_type == "hard")
                        rxLLR_FB = outer_err_seq_est;
                    else
                        rxLLR_FB = 1 - 2*outer_err_seq_est;
                    end
                    inner_err_seq_est = conv_dec(rxLLR_FB,comp_rate,dec_type);
                else
                    inner_err_seq_est = conv_dec(rxLLR_FB,comp_rate,dec_type);
                end
            else
                 % Rate recovery and Decoding
                rxLLR_FB_rr = nrRateRecoverLDPC(rxLLR_FB, K_err, R_err, rv, modulation, nlayers, ncb, Nref);
                [inner_err_seq_est, ~]  = nrldpc_dec(rxLLR_FB_rr, K_err, max_iter, bgn_err);
            end

            % Decompress
            err_deseq_est = arithdeco(double(inner_err_seq_est),counts,length(data))-1;

            % Correct regardless of CRC but store the data (before
            % correction) : except for last round
            data_est_FB_prev = data_est_FB;
            num_err_FB_prev = sum(data ~= double(data_est_FB));

            %% FIX ME : temp 
            data_est_FB_temp = mod(data_est_FB + err_deseq_est,2);
            num_err_FB_temp = sum(data ~= double(data_est_FB_temp));
            if ((i_r < max_rounds - 1) || (i_r == max_rounds -1 && num_err_FB_temp == 0))
                data_est_FB = mod(data_est_FB + err_deseq_est,2);
                num_err_FB = sum(data ~= double(data_est_FB));
            end
            
            % update err count
            num_err_vec(i_r+1) = num_err_FB;
            
            rvSeq_ind = rvSeq_ind + 1;
            
            % Break if crc passes
            if (num_err_FB == 0)
%                 disp(" ************* CORRECTED WITH SPARSE ************* ")
                break;
            end
        else
            % Retransmit
            % generate new tx data
            
            if (code_type == "Conv")         
                [dataIn, rr_len] = conv_enc(data, R);
            else
                % generate new tx data with new rv
                rv = rvSeq(i_r+1);
                
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
            if (channel == "awgn")
                newRxSig = awgn(txSig,SNRdB,'measured');
            elseif (channel == "rayleigh")
                h = sqrt(rand(1)^2 + rand(1)^2);
                newRxSig = awgn(h*txSig,SNRdB,'measured');
            end

            % Symbol demod
            if (modulation == "BPSK")
                if (dec_type == "hard")
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Hard decision';
                    newRxLLR_HARQ = bpskDemodulator(newRxSig);
                else
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Approximate log-likelihood ratio';
                    newRxLLR_HARQ = bpskDemodulator(newRxSig);
                end
            else
                if (dec_type == "hard")
                    newRxLLR_HARQ = qamdemod(newRxSig, M, OutputType='bit', UnitAveragePower=true);
                else
                    newRxLLR_HARQ = qamdemod(newRxSig, M, OutputType='approxllr', UnitAveragePower=true, NoiseVariance=noiseVar);
                end
            end
            
            % Chase Combining
            if (code_type == "LDPC")
                rxLLR_HARQ_buffer = nrRateRecoverLDPC(rxLLR, K, R, rv, modulation, nlayers, ncb, Nref);
                % Rate recovery
                newRxLLR_HARQ = nrRateRecoverLDPC(newRxLLR_HARQ, K, R, rv, modulation, nlayers, ncb, Nref);
            end

            if (dec_type == "hard")
                rxLLR_HARQ_buffer = [rxLLR_HARQ_buffer newRxLLR_HARQ];
                rxLLR_FB_HARQ = round(mean(rxLLR_HARQ_buffer,2));
            else
                rxLLR_HARQ_buffer = [rxLLR_HARQ_buffer newRxLLR_HARQ];
                rxLLR_FB_HARQ = sum(rxLLR_HARQ_buffer,2);
            end
            
            % Rate recovery and Decoding
            if (code_type == "Conv")
                data_est_FB = conv_dec(rxLLR_FB_HARQ, R, dec_type);
            else                
                % Decoding
                bgn = bgn_select(K,R);
                [data_est_FB, crc_chk_FB_HARQ] = nrldpc_dec(rxLLR_FB_HARQ, K, max_iter, bgn);
            end

            % Check for errors : post descrambling
            num_err_FB = sum(mod(data+double(data_est_FB),2));

            % update err count
            num_err_vec(i_r+1) = num_err_FB;
            
            % Break if crc passes
            if (num_err_FB == 0)
                % disp(" ************* CORRECTED WITH HARQ ************* ")
                break;
            end
        end

    end

    % output
    out.Avg_rounds_FB = Avg_rounds_FB;
    out.num_err_FB = num_err_FB;
    out.num_err_vec = num_err_vec;
end