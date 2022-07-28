function out = retransmit_func_FB(SNRdB,modulation,max_iter,rvSeq,nlayers,nPRB,NREPerPRB,N,K,R,data,rxLLR_dsc_rr,ncb,Nref,data_est,err_thr,err_thr_ada_list_est,err_thr_ada_scheme,i_s,max_rounds,counts,num_err,qam_mod,mod_approx,seed)
    
    rng(seed);
    Avg_rounds_FB = 0;
    rxLLR_dsc_rr_HARQ_buffer = rxLLR_dsc_rr;
    rxLLR_FB_dsc_rr_buffer = [];
    data_est_FB = data_est;
    data_est_FB_prev = data_est;
    decision_switch = 0;
    fb_scheme = "HARQ";
    num_err_FB_prev = num_err;
    num_err_vec = zeros(max_rounds,1);
    num_err_vec(1) = num_err;

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
                rxLLR_FB_dsc_rr_buffer = [];
                % Reset rvSeq start
                rvSeq_ind = 1;
            end
            % Use IR for inner scheme as well
            rv = rvSeq(rvSeq_ind); 

            if (inner_scheme == "FB")

                err_seq = arithenco(data_est_err+1,counts);
    
                assert(length(err_seq) < length(data_est_err),'Compression failed in FB round %d',i_r);
                
                targetErrCodeRate = length(err_seq)/N;
    
                % Pick nPRB such that K_err*bits_per_symbol >= length(err_seq)
                nPRB_out = nPRB_select(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate,length(err_seq));
                tbs_err = nPRB_out.tbs_err;
                nPRB = nPRB_out.nPRB;
                
                M = bits_per_symbol(modulation);
                K_err = round(tbs_err/M);
                N_err = nPRB*NREPerPRB;
                R_err = K_err/N_err;
    
                bgn_err = bgn_select(K_err,R_err);
    
                % Padding with zeroes
                nz = K_err - length(err_seq);
                err_seq_n = [err_seq; zeros(nz,1)];
                
                % Encoding and Rate matching
                errDataInSeq = nrldpc_enc(err_seq_n, R_err, modulation, rv, bgn_err, nlayers);

                % Scrambling, TS 38.211 Section 7.3.1.1
                nid = 1; rnti = 1; cwi = 1;
                sc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(errDataInSeq));
                errDataInSeq_sc = xor(errDataInSeq,sc_seq);

                % Symbol Modulation
                if (qam_mod == 1)
                    % Reshape data into binary k-tuples, k = log2(M), and convert to
                    % integers
                    txData = reshape(errDataInSeq_sc,length(dataIn)/M,M);
                    txDataSym = bi2de(txData);
                    txSig_FB = qammod(txDataSym,2^M,'bin','UnitAveragePower',1);
                else
                    txSig_FB = nrSymbolModulate(errDataInSeq_sc,modulation);
                end
                    
                % Pass through AWGN channel
                rxSig_FB = awgn(txSig_FB,SNRdB);
            else
                 % Encoding and Rate matching using new RV
                 errDataInSeq = nrldpc_enc(err_seq_n, R_err, modulation, rv, bgn_err, nlayers);

                 % Scrambling, TS 38.211 Section 7.3.1.1
                 nid = 1; rnti = 1; cwi = 1;
                 sc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(errDataInSeq));
                 errDataInSeq_sc = xor(errDataInSeq,sc_seq);
 
                 % Symbol Modulation
                 if (qam_mod == 1)
                     % Reshape data into binary k-tuples, k = log2(M), and convert to
                     % integers
                     txData = reshape(errDataInSeq_sc,length(dataIn)/M,M);
                     txDataSym = bi2de(txData);
                     txSig_FB = qammod(txDataSym,2^M,'bin','UnitAveragePower',1);
                 else
                     txSig_FB = nrSymbolModulate(errDataInSeq_sc,modulation);
                 end

                % Retransmit through AWGN channel
                rxSig_FB = awgn(txSig_FB,SNRdB);
            end
    
            % Symbol demod
            if (mod_approx == 0)
                noiseVar = 1./(10.^(SNRdB/10));
                rxLLR_FB = nrSymbolDemodulate(rxSig_FB,modulation,noiseVar);
            else
                rxLLR_FB = 1 - 2*double(nrSymbolDemodulate(rxSig_FB,modulation,'DecisionType','hard'));
            end

            % Descrambling, inverse of TS 38.211 Section 7.3.1.1
            opts.MappingType = 'signed';
            opts.OutputDataType = 'double';
            nid = 1; rnti = 1; cwi = 1;
            dsc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(rxLLR_FB),opts);
            rxLLR_FB_dsc = rxLLR_FB .* dsc_seq;

            % Rate recovery
            rxLLR_FB_dsc_rr = nrRateRecoverLDPC(rxLLR_FB_dsc, K_err, R_err, rv, modulation, nlayers, ncb, Nref);

            % Combining
            rxLLR_FB_dsc_rr_buffer = [rxLLR_FB_dsc_rr_buffer rxLLR_FB_dsc_rr];
            rxLLR_FB_dsc_rr = sum(rxLLR_FB_dsc_rr_buffer, 2);

            % Decoding
            [err_seq_est, crc_chk_FB]  = nrldpc_dec(rxLLR_FB_dsc_rr, K_err, max_iter, bgn_err);

            % Decompress
            err_deseq_est = arithdeco(double(err_seq_est),counts,length(data))-1;

            % Correct regardless of CRC but store the data (before
            % correction) : except for last round
            data_est_FB_prev = data_est_FB;
            num_err_FB_prev = sum(data ~= double(data_est_FB));
            if ((i_r < max_rounds - 1) || (i_r == max_rounds -1 && crc_chk_FB == 0))
                data_est_FB = mod(data_est_FB + err_deseq_est,2);
                num_err_FB = sum(data ~= double(data_est_FB));
            end
            
            % update err count
            num_err_vec(i_r+1) = num_err_FB;
            
            rvSeq_ind = rvSeq_ind + 1;
            
            % Break if crc passes
            if (crc_chk_FB == 0)
                break;
            end
        else
            % Retransmit
            % redundancy version for current round
            rv = rvSeq(i_r+1);
            % generate new tx data with new rv
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

            newRxSig = awgn(txSig,SNRdB);

            % Symbol demod
            if (mod_approx == 0)
                noiseVar = 1./(10.^(SNRdB/10));
                newRxLLR_HARQ = nrSymbolDemodulate(newRxSig,modulation,noiseVar);
            else
                newRxLLR_HARQ = 1 - 2*double(nrSymbolDemodulate(newRxSig,modulation,'DecisionType','hard'));
            end
            
            % Descrambling, inverse of TS 38.211 Section 7.3.1.1
            opts.MappingType = 'signed';
            opts.OutputDataType = 'double';
            nid = 1; rnti = 1; cwi = 1;
            dsc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(newRxLLR_HARQ),opts);
            newRxLLR_HARQ_dsc = newRxLLR_HARQ .* dsc_seq;

            % Rate recovery
            newRxLLR_HARQ_dsc_rr = nrRateRecoverLDPC(newRxLLR_HARQ_dsc, K, R, rv, modulation, nlayers, ncb, Nref);
            
            % Combining
            rxLLR_dsc_rr_HARQ_buffer = [rxLLR_dsc_rr_HARQ_buffer newRxLLR_HARQ_dsc_rr];
            rxLLR_FB_HARQ_dsc_rr = sum(rxLLR_dsc_rr_HARQ_buffer,2);
            
            % Decoding
            bgn = bgn_select(K,R);
            [data_est_FB, crc_chk_FB_HARQ] = nrldpc_dec(rxLLR_FB_HARQ_dsc_rr, K, max_iter, bgn);
     
            % Check for errors : post descrambling
            num_err_FB = sum(mod(data+double(data_est_FB),2));

            % update err count
            num_err_vec(i_r+1) = num_err_FB;
            
            % Break if crc passes
            if (crc_chk_FB_HARQ == 0)
                break;
            end
        end

    end

    % output
    out.Avg_rounds_FB = Avg_rounds_FB;
    out.num_err_FB = num_err_FB;
    out.num_err_vec = num_err_vec;
end