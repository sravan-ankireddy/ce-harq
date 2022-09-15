function [output] = conv_dec(llr, rate)

    % quatize the llr to get hard est
    % output = 0;
    llr = llr > 0;
    % if (rate == 1/2)
        trellis = poly2trellis(4,[17 13]);
        tb_len = min(15, round(length(llr)/3));
        decoded = vitdec(llr, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        output = decoded;        
    % end
end