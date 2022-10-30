function rate = converse_rate(eps, N, snr_db, channel, modulation)
% Polyanskiy's achievability bound for VLSF code
    if strcmp(modulation, 'QPSK')
        snr_db = snr_db+ log10(2)*10;
    end
    if strcmp(channel, 'BIAWGN')
        P_by_N = 10.^(snr_db/10);
        capacity = log2(1+P_by_N);
%         disp(capacity);

        logM = (N*capacity + entropy(eps))/(1-eps);

        rate = logM/N;

%         disp(capacity/(1-entropy(eps)));
    end
    disp(capacity/(1-entropy(eps)));
end

function h = entropy(eps) 
    h = -1*(1-eps)*log2(1-eps) - eps*log2(eps);
end
