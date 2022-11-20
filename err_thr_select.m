function err_thr = err_thr_select(data,acomp_table,R,SNRdB,tot_rounds,remRounds,min_bler)

    % Starting with simple scheme : for a given SNR and the round of communication, check
    % if the the error vector is sparse enough to be compressed based on
    % our estimater of err_thr

    % For each SNR, what is the max. rate that can achieve min. ber?

    % Extract the SNR and data
    snr_data = data.snr_data;
    rates_data = data.code_rates;
    
    % find the ind corresponding to current SNR
    [~,snr_ind] = min(abs(snr_data - SNRdB));

    % extract the bler data corresponding to this SNR and the remaining
    % rounds, for all rates
    bler_data = squeeze(data.bler_data(:,remRounds,snr_ind));

    % find the ind corresponding to current rate
    [~,rate_ind] = min(abs(rates_data - R));

    % bler corresponding to max_rounds of HARQ with base rate
    bler_data_R = squeeze(data.bler_data(rate_ind,tot_rounds,snr_ind));
    % disp ("base bler : ");
    % disp(bler_data');
    % disp("new bler : rem rounds ");
    % disp(remRounds);
    % disp(bler_data_R');
    % find all the rates with bler < min_bler
    rates_ind_list = find(bler_data <= bler_data_R);
    % rates_ind_list = find(bler_data <= min_bler);
    if (~isempty(rates_ind_list))
    
    % bler_data
    % bler_data_R
    % min_bler
    % rates_ind_list
    % if (min(bler_data) > 1e-2)
    % if (bler_data_R > 1e-2)
    %     err_thr = 0.05;
    % % else
    % elseif (length(rates_ind_list) >= 1)
        % find the max rate
        max_rate = max(rates_data(rates_ind_list));

        % use this rate to calculate the rate of compression needed
        Rc = max_rate/0.75;

        % use the compression LUT to pick the sparsity 
        rc_data = acomp_table.Rc_vec;
        sp_data = acomp_table.sp_vec;

        % get list of comp rates < Rc
        rc_data_feasible = rc_data(rc_data < Rc);

        [~, rc_ind] = min(abs(rc_data_feasible - Rc));
        
        err_thr = sp_data(rc_ind);

        % if (SNRdB == 1.6)
        %     disp("R")
        %     disp(R)
        %     disp("rem rounds")
        %     disp(remRounds);
        %     disp("err thr")
        %     disp(err_thr);
        %     pause(1)
        % end
    else
        err_thr = 0;
    end

    
end