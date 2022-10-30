function out = retransmit_func_FB_vit(channel,SNRdB,modulation,N,K,R,data,rxLLR,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,comm_mod,mod_approx,seed)
    
    rng(seed);
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

    rep_coding = 0;
    tc_rp = 0.2;
    rc = 0;

    comp_rates = 1./(2:12);
    minR = comp_rates(end);

    for i_r = 1:max_rounds-1

        num_err_FB = sum(data ~= double(data_est_FB));
        err_per = num_err_FB/K;

        % Adaptive err thr
        if (err_thr_ada_scheme == "est")
            err_thr = err_thr_ada_list_est((max_rounds - i_r),i_s);
        end

        % Once the error becomes sparse enough, stay on FB scheme
        % Within the FB scheme, you can choose FB/HARQ based on error
        % vector
        if (decision_switch == 0)
            if (err_per <= err_thr && err_thr > 0)
                % disp("found sparse err ")
                % Also check if the error is compressible
                data_est_err_temp = mod(data+double(data_est_FB),2);
                err_seq_temp = arithenco(data_est_err_temp+1,counts);
                if (length(err_seq_temp) < length(data_est_err_temp))
                    fb_scheme = "FB";
                    decision_switch = 1;
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
                % Empty the chase buffer if prev round reduced the errors

                rxLLR_FB_buffer = [];
                % Reset rvSeq start
                rvSeq_ind = 1;
            end

            if (inner_scheme == "FB")
                % disp("inner feedback scheme ")

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
                    min_ind = length(comp_rates) - 1;
                end

                comp_str = 'conv' + string(min_ind+1);
                
                % encode to length <= K
                comp_code = channel_encode(err_seq, comp_str);

%                 % pad with zeros
%                 nnz = K - length(comp_code);
%                 assert(nnz >= 0,'Comp code failed in FB round %d',i_r);
% 
%                 comp_code = [comp_code; zeros(nnz,1)];

                errDataInSeq = comp_code;
                
%                 % use rep coding to go beyond least code rate
%                 if (rep_coding == 1 && targetErrCodeRate < tc_rp)
%                     rc = floor(tc_rp*N/(length(err_seq)+16));
%                     err_seq_vec = [];
%                     for ie = 1:rc
%                         err_seq_vec = [err_seq_vec; err_seq];
%                     end
%                     err_seq = err_seq_vec;
%                     targetErrCodeRate = (length(err_seq)+16)/N;
%                 end

                % Symbol Modulation
                bpskModulator = comm.BPSKModulator;
                bpskModulator.PhaseOffset = pi/4;
                txSig_FB = bpskModulator(errDataInSeq);
                    
                % Pass through channel
                if (channel == "awgn")
                    rxSig_FB = awgn(txSig_FB,SNRdB);
                elseif (channel == "rayleigh")
                    h = sqrt(rand(1)^2 + rand(1)^2);
                    rxSig_FB = awgn(h*txSig_FB,SNRdB);
                end
            else
                % Encoding and Rate matching using new RV
%                 [errDataInSeq, err_rr_len] = conv_enc(inner_code_int, R);
 
                % Symbol Modulation
                bpskModulator = comm.BPSKModulator;
                bpskModulator.PhaseOffset = pi/4;
                txSig_FB = bpskModulator(errDataInSeq);

                % Pass through channel
                if (channel == "awgn")
                    rxSig_FB = awgn(txSig_FB,SNRdB);
                elseif (channel == "rayleigh")
                    h = sqrt(rand(1)^2 + rand(1)^2);
                    rxSig_FB = awgn(h*txSig_FB,SNRdB);
                end
            end
    
            % Symbol demod
            if (comm_mod == 1)
                if (mod_approx == 0)
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Approximate log-likelihood ratio';
                    rxLLR_FB = -1*bpskDemodulator(rxSig_FB);
                else
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Hard decision';
                    rxLLR_FB = -1 + 2*bpskDemodulator(rxSig_FB);
                end
            else
                if (mod_approx == 0)
                    noiseVar = 1./(10.^(SNRdB/10));
                    rxLLR_FB = lteSymbolDemodulate(rxSig_FB,modulation,noiseVar);
                else
                    rxLLR_FB = 1 - 2*double(lteSymbolDemodulate(rxSig_FB,modulation,'DecisionType','hard'));
                end
            end
            
%             % Rep code combining
%             if (rc > 1 && rep_coding == 1)
%                 temp = reshape(rxLLR_FB, length(rxLLR_FB)/rc,rc);
%                 rxLLR_FB = sum(temp,2);
%             end

            % Combining
            rxLLR_FB_buffer = [rxLLR_FB_buffer rxLLR_FB];
            rxLLR_FB = sum(rxLLR_FB_buffer, 2);
            
%             % Outer code : decode outer code
%             outer_err_seq_est = double(conv_dec(rxLLR_FB, err_rr_len));
% 
%             % Deinterleaver            
%             outer_err_seq_est = randdeintrlv(outer_err_seq_est,0);
%     
%             % Inner code : decoder inner coder using llr from hard info 
% %             outer_llr = 2*outer_err_seq_est - 1;
%             % remove the padded zeros
%             outer_err_seq_est = outer_err_seq_est(1:end-nnz);
            outer_err_seq_est = rxLLR_FB > 0;
            inner_err_seq_est = channel_decode(outer_err_seq_est,comp_str);

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
                disp(" ************* CORRECTED WITH SPARSE ************* ")
                break;
            end
        else
            % Retransmit
            % generate new tx data         
            [dataIn, rr_len] = conv_enc(data, R);

            % Symbol Modulation
            if (comm_mod == 1)
                bpskModulator = comm.BPSKModulator;
                bpskModulator.PhaseOffset = pi/4;
                txSig = bpskModulator(dataIn);
            else
                txSig = lteSymbolModulate(dataIn,modulation);
            end

            % Pass through channel
            if (channel == "awgn")
                newRxSig = awgn(txSig,SNRdB);
            elseif (channel == "rayleigh")
                h = sqrt(rand(1)^2 + rand(1)^2);
                newRxSig = awgn(h*txSig,SNRdB);
            end

            % Symbol demod
            if (comm_mod == 1)
                if (mod_approx == 0)
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Approximate log-likelihood ratio';
                    newRxLLR_HARQ = -1*bpskDemodulator(newRxSig);
                else
                    bpskDemodulator = comm.BPSKDemodulator; 
                    bpskDemodulator.PhaseOffset = pi/4; 
                    bpskDemodulator.DecisionMethod = 'Hard decision';
                    newRxLLR_HARQ = -1 + 2*bpskDemodulator(newRxSig);
                end
            else
                if (mod_approx == 0)
                    noiseVar = 1./(10.^(SNRdB/10));
                    newRxLLR_HARQ = lteSymbolDemodulate(newRxSig,modulation,noiseVar);
                else
                    newRxLLR_HARQ = 1 - 2*double(lteSymbolDemodulate(newRxSig,modulation,'DecisionType','hard'));
                end
            end
            
            % Combining
            rxLLR_HARQ_buffer = [rxLLR_HARQ_buffer newRxLLR_HARQ];
            rxLLR_FB_HARQ = sum(rxLLR_HARQ_buffer,2);
            
            % Rate recovery and Decoding
            data_est_FB = conv_dec(rxLLR_FB_HARQ, rr_len);
     
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