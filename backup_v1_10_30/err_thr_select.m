function err_thr = err_thr_select(data,acomp_table,R,SNRdB,remRounds,min_bler)

    % Starting with simple scheme : for a given SNR and the round of communication, check
    % if the the error vector is sparse enough to be compressed based on
    % our estimater of err_thr

    % For each SNR, what is the max. rate that can achieve min. ber?

    % Extract the SNR and data
    snr_data = data.snr_data;
    rates_data = data.actualCodeRate_list;
    
    % find the ind corresponding to current SNR
    [~,snr_ind] = min(abs(snr_data - SNRdB));

    % extract the bler data corresponding to this SNR and the remaining
    % rounds, for all rates
    bler_data = squeeze(data.err_data_all_rates(2,:,remRounds,snr_ind));

    % find all the rates with bler < min_bler
    rates_ind_list = find(bler_data < min_bler);

    if (length(rates_ind_list) < 1)
        rates_ind_list = 1;
    end

    % find the max rate
    max_rate = rates_data(max(rates_ind_list));

    % use this rate to calculate the rate of compression needed
    Rc = max_rate/R;

    % use the compression LUT to pick the sparsity 
    rc_data = acomp_table.Rc_vec;
    sp_data = acomp_table.sp_vec;
    [~, rc_ind] = min(abs(rc_data - Rc));
    
    err_thr = sp_data(rc_ind);
end