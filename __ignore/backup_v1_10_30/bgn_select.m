function bgn = bgn_select(K,R)
    if (K <= 3824 && R <= 0.67)
        bgn = 2;
    elseif (K <= 292)
        bgn = 2;
    elseif (R <= 0.25)
        bgn = 2;
    else
        bgn = 1;
    end
end


