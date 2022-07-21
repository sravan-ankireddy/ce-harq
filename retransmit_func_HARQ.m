function out = retransmit_func_HARQ(SNRdB,modulation,max_iter,rv,nlayers,K,R,data,txSig,rxLLR,err_thr,num_err,max_rounds,mod_approx,seed)
    
    rng(seed);

    % Init
    Avg_rounds_HARQ = 0;

    % Avg rounds post sparsity
    Avg_rounds_HARQ_s = 0;
    
    num_err_HARQ = inf;
    num_err_vec = zeros(max_rounds,1);
    num_err_vec(1) = num_err;
    
    %% Always run HARQ for baseline
    rxLLR_HARQ_chase = rxLLR;
    for i_r = 1:max_rounds-1
        
        err_per = num_err_HARQ/K;
        if (err_per < err_thr)
            Avg_rounds_HARQ_s = Avg_rounds_HARQ_s + 1;
        end

        % Update the counter
        Avg_rounds_HARQ = Avg_rounds_HARQ + 1; 
        newRxSig = awgn(txSig,SNRdB);

        % Symbol demod
        if (mod_approx == 0)
            noiseVar = 1./(10.^(SNRdB/10));
            newRxLLR = nrSymbolDemodulate(newRxSig,modulation,noiseVar);
        else
            newRxLLR = 1 - 2*double(nrSymbolDemodulate(newRxSig,modulation,'DecisionType','hard'));
        end

        % Chase combining
        rxLLR_HARQ_chase = [rxLLR_HARQ_chase newRxLLR];
        rxLLR_HARQ = sum(rxLLR_HARQ_chase,2);

        % Descrambling, inverse of TS 38.211 Section 7.3.1.1
        opts.MappingType = 'signed';
        opts.OutputDataType = 'double';
        nid = 1; rnti = 1; cwi = 1;
        dsc_seq = nrPDSCHPRBS(nid,rnti,cwi-1,length(rxLLR_HARQ),opts);
        rxLLR_HARQ_dsc = rxLLR_HARQ .* dsc_seq;
                
        % Decoding and Rate recovery
        bgn = bgn_select(K,R);
        [data_est_HARQ, crc_chk_HARQ]  = nrldpc_dec(rxLLR_HARQ_dsc, R, modulation, K, max_iter, rv, nlayers, bgn);

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