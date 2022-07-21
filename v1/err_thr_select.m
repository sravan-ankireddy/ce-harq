function err_thr = err_thr_select(N,R,SNR, min_ber)
    % For each SNR, what is the max. rate required to get a BER of 10^-5?
    filename = sprintf('results_rate_curves/BER_LDPC_%d.mat',N);
    load(filename,'data');
    
    snr_vec = round(data(1,:),4);
    snr_ind = find(snr_vec==SNR);

    rate_vec = round(data(2:end,1),4);

    % BER vec for that rate
    BER_vec = data(2:end,snr_ind(1));
    
    max_rate_ind = -1;
    for i = 1:length(BER_vec)
        if BER_vec(i) < min_ber
            max_rate_ind = i;
        end
    end    
    
    if (max_rate_ind > 0)
        max_rate = rate_vec(max_rate_ind);
        % Max compression rate
        Rc = max_rate/R;
    else
        Rc = 0.1;
    end

    % Use avg. compression rate to calculate the max allowed errors
%     [Avg_comp, sp_vec] = acomp_table(N);
    out = load('results_rate_curves/acomp_960.mat');
    Avg_comp = out.Rc_vec;
    sp_vec = out.sp_vec;
    
    max_comp_ind = -1;
    for j = 1:length(sp_vec)
        if (Avg_comp(j) < Rc)
            max_comp_ind = j;
        end
    end

    err_thr = sp_vec(max_comp_ind);
end