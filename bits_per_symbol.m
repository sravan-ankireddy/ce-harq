function bps = bits_per_symbol(m)
    
    if (m == "BPSK" || m == "pi/2-BPSK")
        bps = 1;
    elseif(m == "QPSK")
        bps = 2;
    elseif(m == "16QAM")
        bps = 4;
    elseif(m == "64QAM")
        bps = 6;
    elseif(m == "256QAM")
        bps = 8;
    else
        bps = -1;
    end
end