function out = retransmit_func_HARQ(SNRdB,modulation,max_iter,rvSeq,nlayers,K,R,data,rxLLR_dsc_rr,ncb,Nref,err_thr,num_err,max_rounds,qam_mod,mod_approx,seed)
    
    rng(seed);

    % Init
    Avg_rounds_HARQ = 0;

    % Avg rounds post sparsity
    Avg_rounds_HARQ_s = 0;
    
    num_err_HARQ = inf;
    num_err_vec = zeros(max_rounds,1);
    num_err_vec(1) = num_err;
    
    %% Buffer to store rxLLR after descarmbling and rate recovery
    rxLLR_HARQ_dsc_rr_buffer = rxLLR_dsc_rr;
    for i_r = 1:max_rounds-1
        % redundancy version for current round
        rv = rvSeq(i_r+1);
        
        err_per = num_err_HARQ/K;
        if (err_per < err_thr)
            Avg_rounds_HARQ_s = Avg_rounds_HARQ_s + 1;
        end

        % Update the counter
        Avg_rounds_HARQ = Avg_rounds_HARQ + 1; 
        
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
            newRxLLR = nrSymbolDemodulate(newRxSig,modulation,noiseVar);
        else
            newRxLLR = 1 - 2*double(nrSymbolDemodulate(newRxSig,modulation,'DecisionType','hard'));
        end

        % Descrambling, inverse of TS 38.211 Section 7.3.1.1
        opts.MappingType = 'signed';
        opts.OutputDataType = 'double';
        nid = 1; rnti = 1; cwi = 1;
        dsc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(newRxLLR),opts);
        newRxLLR_dsc = newRxLLR.* dsc_seq;
        
        % Rate recovery
        newRxLLR_dsc_rr = nrRateRecoverLDPC(newRxLLR_dsc, K, R, rv, modulation, nlayers, ncb, Nref);

        % Combining
        rxLLR_HARQ_dsc_rr_buffer = [rxLLR_HARQ_dsc_rr_buffer newRxLLR_dsc_rr];
        rxLLR_HARQ_dsc_rr = sum(rxLLR_HARQ_dsc_rr_buffer,2);
                
        % Decoding
        bgn = bgn_select(K,R);
        [data_est_HARQ, crc_chk_HARQ]  = nrldpc_dec(rxLLR_HARQ_dsc_rr, K, max_iter, bgn);

        % Check for errors
        num_err_HARQ = sum(data ~= double(data_est_HARQ));

        num_err_vec(i_r+1) = num_err_HARQ;

        if (crc_chk_HARQ == 0)
            break;
        end
    end

    % output
    out.Avg_rounds_HARQ = Avg_rounds_HARQ;
    out.Avg_rounds_HARQ_s = Avg_rounds_HARQ_s;
    out.num_err_HARQ = num_err_HARQ;
    out.num_err_vec = num_err_vec;
end