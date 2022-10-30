function rate = achievable_rate(eps, N, snr_db, channel, modulation)
% Polyanskiy's achievability bound for VLSF codes
    if strcmp(modulation, 'QPSK')
        snr_db = snr_db + 3;
    end
    if strcmp(channel, 'BIAWGN')
        P_by_N = 10.^(snr_db/10);
        capacity = 0.5*log2(1+P_by_N);

        logM = N*capacity + log2(eps) - 1;
        rate = logM/N;
    end
end