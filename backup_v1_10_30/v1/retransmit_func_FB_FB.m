function out = retransmit_func_FB_FB(SNRdB,modulation,max_iter,rv,nlayers,nPRB,NREPerPRB,M,N,K,R,data,txSig,rxSig,data_est,err_thr,max_rounds,counts,num_err)
    
    Avg_rounds_FB = 0;
    rxSig_FB_HARQ_chase = rxSig;
    rxSig_FB_chase = [];
    data_est_FB = data_est;
    data_est_FB_prev = data_est;
    decision_switch = 0;
    fb_scheme = "HARQ";
    num_err_FB_prev = num_err;
    
    for i_r = 1:max_rounds-1
        num_err_FB = sum(mod(data+double(data_est_FB),2));
        err_per = num_err_FB/K;

        % Once the error becomes sparse enough, stay on FB scheme
        % Within the FB scheme, you can choose FB/HARQ based on error
        % vector
        if (decision_switch == 0)
            if (err_per <= err_thr)
                % Also check if the error is compressible
                data_est_err_temp = mod(data+double(data_est_FB),2);
                err_seq_temp = arithenco(data_est_err_temp+1,counts);
                if (length(err_seq_temp) < length(data_est_err_temp))
                    fb_scheme = "FB";
                    decision_switch = 1;
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
            else
                inner_scheme = "FB";
                % Empty the chase buffer if prev round reduced the errors
                rxSig_FB_chase = [];
            end

            if (inner_scheme == "FB")

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
                
                % Store for chase combining
                rxSig_FB_chase = [rxSig_FB_chase rxSig_FB];

            else
                % Retransmit through AWGN channel
                rxSig_FB = awgn(txSig_FB,SNRdB);

                % Store for chase combining
                rxSig_FB_chase = [rxSig_FB_chase rxSig_FB];
                rxSig_FB = mean(rxSig_FB_chase, 2);
            end
    
            % QAM Demod
            rxLLR_FB = qamdemod(rxSig_FB,M,'OutputType','LLR');
        
            % Decoding and Rate recovery
            [err_seq_est, crc_chk_FB]  = nrldpc_dec(rxLLR_FB, R_err, modulation, K_err, max_iter, rv, nlayers, bgn_err);

            % Decompress
            err_deseq_est = arithdeco(double(err_seq_est),counts,length(data))-1;

            % Correct regardless of CRC but store the data (before
            % correction)
            data_est_FB_prev = data_est_FB;
            num_err_FB_prev = sum(mod(data+double(data_est_FB),2));

            data_est_FB = mod(data_est_FB + err_deseq_est,2);
            num_err_FB = sum(mod(data+double(data_est_FB),2));

            % Break if crc passes
            if (crc_chk_FB == 0)
                break;
            end

        else
            % Retransmit
            newRxSig = awgn(txSig,SNRdB);

            % Chase combining
            rxSig_FB_HARQ_chase = [rxSig_FB_HARQ_chase newRxSig];
            rxSig_FB_HARQ = mean(rxSig_FB_HARQ_chase,2);

            % QAM Demod
            rxLLR_FB_HARQ = qamdemod(rxSig_FB_HARQ,M,'OutputType','LLR');
                    
            % Decoding and Rate recovery
            bgn = bgn_select(K,R);
            [data_est_FB, crc_chk_FB_HARQ] = nrldpc_dec(rxLLR_FB_HARQ, R, modulation, K, max_iter, rv, nlayers, bgn);
    
            % Check for errors
            num_err_FB = sum(mod(data+double(data_est_FB),2));
            
            % Break if crc passes
            if (crc_chk_FB_HARQ == 0)
                break;
            end

        end

    end

    % output
    out.Avg_rounds_FB_FB = Avg_rounds_FB;
    out.num_err_FB_FB = num_err_FB;
end