function out = retransmit_func_ARQ(SNRdB,modulation,max_iter,rv,nlayers,nPRB,NREPerPRB,M,N,K,R,data,txSig,rxSig,data_est,err_thr,max_rounds,counts,num_err)
    
    % Init
    Avg_rounds_ARQ = 0;


    %% Always run ARQ for baseline
    for i_r = 1:max_rounds-1

        % Update the counter
        Avg_rounds_ARQ = Avg_rounds_ARQ + 1; 
        newRxSig = awgn(txSig,SNRdB);
        rxSig_ARQ = newRxSig;

        % QAM Demod
        rxLLR_ARQ = qamdemod(rxSig_ARQ,M,'OutputType','LLR');
                
        % Decoding and Rate recovery
        bgn = bgn_select(K,R);
        [data_est_ARQ, ~]  = nrldpc_dec(rxLLR_ARQ, R, modulation, K, max_iter, rv, nlayers, bgn);

        % Check for errors
        num_err_ARQ = sum(mod(data+double(data_est_ARQ),2));
        
        if (num_err_ARQ == 0)
            break;
        end
    end

    % output
    out.Avg_rounds_ARQ = Avg_rounds_ARQ;

    out.num_err_ARQ = num_err_ARQ;
end