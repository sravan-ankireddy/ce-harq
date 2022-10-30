function out = nPRB_select(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate,err_seq_len)
    while 1 
        if (nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate) < err_seq_len)
            nPRB = nPRB + 1;
        else
            break;
        end
    end
    out = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetErrCodeRate);
end