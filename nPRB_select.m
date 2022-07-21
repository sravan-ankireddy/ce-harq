function out = nPRB_select(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate,err_seq_len)
    while 1 
        if (nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate)/bits_per_symbol(modulation) < err_seq_len)
            nPRB = nPRB + 1;
        else
            break;
        end
    end
    out.tbs_err = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate);
    out.nPRB = nPRB;
end