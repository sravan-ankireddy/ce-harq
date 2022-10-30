function [output, rr_len] = conv_enc(data_in, rate)
    
    output = 0;
    if (rate == 5/6)
        trellis = poly2trellis(4,[17 13]);
        puncpat = [1;1;0;1;0];
        % FIX ME
        if (round(length(data_in)/rate) <= 400)
            termination = [0, 0, 0, 0, 0, 0, 0];
        else
            termination = [0, 0, 0, 0];
        end
        temp = [data_in; termination'];
        output = convenc([data_in; termination'], trellis, puncpat);
    elseif (rate == 3/4)
        trellis = poly2trellis(4,[17 13]);
        termination = [0, 0, 0];
        puncpat = [1;1;0];
        output = convenc([data_in; termination'], trellis, puncpat);
    elseif (rate == 1/2)
        trellis = poly2trellis(4,[17 13]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/3)
        trellis = poly2trellis(4,[17 15 13]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/4)
        trellis = poly2trellis(4,[17 15 13 11]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/5)
        trellis = poly2trellis(4,[17 15 13 11 7]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/6)
        trellis = poly2trellis(4,[17 15 13 11 7 5]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/7)
        trellis = poly2trellis(4,[17 15 13 11 7 5 3]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/8)
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/9)
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/10)
        trellis = poly2trellis(4,[17 15 13 13 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/11)
        trellis = poly2trellis(4,[17 15 13 13 11 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    elseif (rate == 1/12)
        trellis = poly2trellis(4,[17 15 15 13 13 11 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        output = convenc([data_in; termination'], trellis);
    end
    rr_len = length(output);
end