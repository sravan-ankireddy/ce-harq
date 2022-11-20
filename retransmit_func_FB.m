function out = retransmit_func_FB(channel,SNRdB,modulation,N,K,R,dec_type,data,rxLLR,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,comm_mod,mod_approx,seed)
    
    rng(seed);
    M = bits_per_symbol(modulation);
    Avg_rounds_FB = 0;
    rxLLR_HARQ_buffer = rxLLR;
    rxLLR_FB_buffer = [];
    data_est_FB = data_est;
    data_est_FB_prev = data_est;
    decision_switch = 0;
    fb_scheme = "HARQ";
    num_err_FB_prev = num_err;
    num_err_vec = zeros(max_rounds,1);
    num_err_vec(1) = num_err;

    noiseVar = 1./(10.^(SNRdB/10));
    rep_coding = 0;

    comp_rates = 1./(2:12);
    minR = comp_rates(end);

    min_bler = 1e-4;
    harq_data = load('bler_data/awgn/400/hard/harq/BPSK/10000/harq_data_Conv_400_rate_0.833_rate_0.083_max_rounds_4.mat');
    acomp_table = load('lut_data/acomp_400_ns_100000.mat');
    prev_reset_round = max_rounds;
    for i_r = 1:max_rounds-1

        num_err_FB = sum(data ~= double(data_est_FB));
        err_per = num_err_FB/K;

        % Adaptive err thr
        if (err_thr_ada_scheme == "est")
            % err_thr = err_thr_ada_list_est((max_rounds - i_r),i_s);
            if (fb_scheme == "FB")
                err_thr = err_thr_select(harq_data,acomp_table,targetErrCodeRate,SNRdB,max_rounds - prev_reset_round+1,max_rounds - i_r,min_bler);
            else
                err_thr = err_thr_select(harq_data,acomp_table,R,SNRdB,max_rounds,max_rounds - i_r,min_bler);
            end
            err_thr_ada_list_est((max_rounds - i_r),i_s) = err_thr;
        end
        % err_tr = 0.05;
        % Once the error becomes sparse enough, stay on FB scheme
        % Within the FB scheme, you can choose FB/HARQ based on error
        % vector
        if (decision_switch == 0)
            if (err_per <= err_thr && err_thr > 0)
                % Also check if the error is compressible
                data_est_err_temp = mod(data+double(data_est_FB),2);

                % iteratively update counts
                p_cur = max(1,round(100 * sum(data_est_err_temp)/length(data_est_err_temp)));
                counts_iter = [100-p_cur p_cur];
                err_seq_temp = arithenco(data_est_err_temp+1,counts_iter);

                % err_seq_temp = arithenco(data_est_err_temp+1,counts);
                if (length(err_seq_temp) < length(data_est_err_temp))
                    fb_scheme = "FB";
                    decision_switch = 1;
                    counts = counts_iter;
                else
                    fb_scheme = "HARQ";
                end
            else
                fb_scheme = "HARQ";
            end
        end

        % Update the counter
        Avg_rounds_FB = Avg_rounds_FB + 1; 
        rvSeq_ind = 1;
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
                prev_reset_round = i_r;
                % Empty the chase buffer if prev round reduced the errors

                rxLLR_FB_buffer = [];
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

                % find the smallest rate for compression
                targetErrCodeRate = (length(err_seq)+3)/N;

                % find the rate for inner code
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
                end
                
                % encode to length <= N
                comp_code = conv_enc(err_seq, comp_rate);

                errDataInSeq = comp_code;

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
                % Encoding and Rate matching using new RV
%                 [errDataInSeq, err_rr_len] = conv_enc(inner_code_int, R);
 
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

            
            % outer_err_seq_est = rxLLR_FB > 0;
            inner_err_seq_est = conv_dec(rxLLR_FB,comp_rate,dec_type);

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
            [dataIn, rr_len] = conv_enc(data, R);

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
            if (dec_type == "hard")
                rxLLR_HARQ_buffer = [rxLLR_HARQ_buffer newRxLLR_HARQ];
                rxLLR_FB_HARQ = round(mean(rxLLR_HARQ_buffer,2));
            else
                rxLLR_HARQ_buffer = [rxLLR_HARQ_buffer newRxLLR_HARQ];
                rxLLR_FB_HARQ = sum(rxLLR_HARQ_buffer,2);
            end
            
            % Rate recovery and Decoding
            data_est_FB = conv_dec(rxLLR_FB_HARQ, R, dec_type);
     
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