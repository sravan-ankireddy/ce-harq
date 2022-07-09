function out = retransmit_func_FB_HARQ(SNRdB,modulation,max_iter,rv,nlayers,nPRB,NREPerPRB,M,N,K,R,data,txSig,rxSig,data_est,err_thr,max_rounds,counts,num_err)
    
    % CRC pass == successful decoding
    Avg_rounds_FB = 0;
    % buffer for storing original data
    rxSig_buff = rxSig;
    % buffer for storing error vector data
    rxSig_FB_buffer = [];
    data_est_FB = data_est;
    decision_switch = 0;
    fb_scheme = "HARQ";

    for i_r = 1:max_rounds-1
        num_err_FB = sum(mod(data+double(data_est_FB),2));
        err_per = num_err_FB/K;

        % Once the error becomes sparse enough, stay on FB scheme and
        % perform HARQ till error vector is decoded correctly
        if (decision_switch == 0)
            if (err_per <= err_thr)
                % Check if error is compressible
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

            % Start storing the rxSig for HARQ
            rxSig_FB_buffer = [rxSig_FB_buffer rxSig_FB];
            rxSig_FB = mean(rxSig_FB_buffer,2);
    
            % QAM Demod
            rxLLR_FB = qamdemod(rxSig_FB,M,'OutputType','LLR');
        
            % Decoding and Rate recovery
            [err_seq_est, crc_chk_FB]  = nrldpc_dec(rxLLR_FB, R_err, modulation, K_err, max_iter, rv, nlayers, bgn_err);

            % Decompress
            err_deseq_est = arithdeco(double(err_seq_est),counts,length(data))-1;
            
            % Correct the error only if CRC passes else round
            % wasted`
            if (crc_chk_FB == 0)
                data_est_FB = mod(data_est_FB + err_deseq_est,2);
                num_err_FB = sum(data ~= double(data_est_FB));
                break;
            end

        else
            % Retransmit
            newRxSig = awgn(txSig,SNRdB);

            % Chase combining
            rxSig_buff = [rxSig_buff newRxSig];
            rxSig_FB = mean(rxSig_buff,2);

            % QAM Demod
            rxLLR_FB = qamdemod(rxSig_FB,M,'OutputType','LLR');
                    
            % Decoding and Rate recovery
            bgn = bgn_select(K,R);
            [data_est_FB, crc_chk_FB] = nrldpc_dec(rxLLR_FB, R, modulation, K, max_iter, rv, nlayers, bgn);
    
            % Check for errors
            num_err_FB = sum(data ~= double(data_est_FB));

            if (crc_chk_FB == 0)
                break;
            end
        end

    end

    % output
    out.Avg_rounds_FB_HARQ = Avg_rounds_FB;
    out.num_err_FB_HARQ = num_err_FB;
end