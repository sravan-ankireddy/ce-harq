function [output, rr_len] = conv_enc(data_in, rate)
    
    % output = 0;
    % rr_len = 0;
    % if (rate == 1/2)
        trellis = poly2trellis(4,[17 13]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
        rr_len = length(output);        
    % end
end