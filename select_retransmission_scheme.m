function scheme = select_retransmission_scheme(data,R,Re,SNRdB,numR_HARQ,numR_FB)

    % simple scheme : At every check point, compare Prob of error for
    % numR_harq rounds of HARQ ar rate R vs numR_harq rounds of FB at rate
    % Re

    % Extract snr and rate info
    snr_data = data.snr_data;
    rate_data = data.actualCodeRate_list;

    [~,snr_ind] = min(abs(snr_data - SNRdB));
    
    [~,rate_ind_HARQ] = min(abs(rate_data - R));

    % get the corresponding bler for HARQ
    bler_HARQ = data.err_data_all_rates(2,rate_ind_HARQ,numR_HARQ,snr_ind);

    [~,rate_ind_FB] = min(abs(rate_data - Re));

    % get the corresponding bler for FB(HARQ)
    bler_FB = data.err_data_all_rates(2,rate_ind_FB,numR_FB,snr_ind);

    % pick the scheme with lesser expected error
    if (bler_FB <= bler_HARQ)
        scheme = "FB";
    else
        scheme = "HARQ";
    end
end