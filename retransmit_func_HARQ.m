function out = retransmit_func_HARQ(SNRdB,modulation,max_iter,rv,nlayers,M,K,R,data,txSig,rxSig,max_rounds)
    
    % Init
    Avg_rounds_HARQ = 0;

    %% Always run HARQ for baseline
    rxSig_FB_HARQ_chase = rxSig;
    for i_r = 1:max_rounds-1

        % Update the counter
        Avg_rounds_HARQ = Avg_rounds_HARQ + 1; 
        newRxSig = awgn(txSig,SNRdB);

        % Chase combining
        rxSig_FB_HARQ_chase = [rxSig_FB_HARQ_chase newRxSig];
        rxSig_HARQ = mean(rxSig_FB_HARQ_chase,2);

        % QAM Demod
        rxLLR_HARQ = qamdemod(rxSig_HARQ,M,'OutputType','LLR');
                
        % Decoding and Rate recovery
        bgn = bgn_select(K,R);
        [data_est_HARQ, ~]  = nrldpc_dec(rxLLR_HARQ, R, modulation, K, max_iter, rv, nlayers, bgn);

        % Check for errors
        num_err_HARQ = sum(mod(data+double(data_est_HARQ),2));
        
        if (num_err_HARQ == 0)
            break;
        end
    end

    % output
    out.Avg_rounds_HARQ = Avg_rounds_HARQ;

    out.num_err_HARQ = num_err_HARQ;
end