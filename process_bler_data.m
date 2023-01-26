function [bler_out, ar_out, snr_out] = process_bler_data(data)

    snr_data = data.snr_data;
    err_thr_list = data.err_thr_grid;

    % get the bler_data for last round
    bler_data = squeeze(data.bler_data(:,end,:));
    
    bler_out.err_thr_opt = zeros(length(snr_data),1);
    bler_out.bler_opt = zeros(length(snr_data),1);

    % for each SNR, find the err_thr where min. bler os observed
    for is = 1:length(snr_data)
        [min_val, min_ind] = min(bler_data(:,is));

        % Pick the largest err thr if more than one err_thr gives same bler, including 0
        min_ind_list = find(bler_data(:,is) == min_val);
        if (length(min_ind_list) > 1)
            min_ind = max(min_ind_list);
        end
        bler_out.err_thr_opt(is) = err_thr_list(min_ind);
        bler_out.bler_opt(is) = bler_data(min_ind,is);
    end

    % get the ar_data
    ar_data = data.ar_data;
    
    ar_out.err_thr_opt = zeros(length(snr_data),1);
    ar_out.ar_opt = zeros(length(snr_data),1);

    % for each SNR, find the err_thr where min. ar observed
    for is = 1:length(snr_data)
        [min_val, min_ind] = min(ar_data(:,is));

        % Pick the highest err thr if more than one err_thr gives same ar
        min_ind_list = find(ar_data(:,is) == min_val);
        if (length(min_ind_list) > 1)
            min_ind = max(min_ind_list);
        end

        ar_out.err_thr_opt(is) = err_thr_list(min_ind);
        ar_out.ar_opt(is) = ar_data(min_ind,is);
    end
    snr_out = snr_data;
end