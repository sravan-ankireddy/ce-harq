function tx_code = channel_encode(msg, code_mode, interleaver)
    if nargin > 2
        intrlvrIndices = interleaver;
    else
        intrlvrIndices = 0;
    end

    if strcmp(code_mode, 'uncoded')
        tx_code = msg;
    elseif strcmp(code_mode, 'turbo')
        trellis = poly2trellis(4,[13 15],13);
        turboenc = comm.TurboEncoder(trellis,intrlvrIndices);
        tx_code = turboenc(msg')';
    elseif strcmp(code_mode, 'conv2')
        trellis = poly2trellis(4,[17 13]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv3')
        trellis = poly2trellis(4,[17 15 13]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv4')
        trellis = poly2trellis(4,[17 15 13 11]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv5')
        trellis = poly2trellis(4,[17 15 13 11 7]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv6')
        trellis = poly2trellis(4,[17 15 13 11 7 5]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv7')
        trellis = poly2trellis(4,[17 15 13 11 7 5 3]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv8')
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv9')
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv10')
        trellis = poly2trellis(4,[17 15 13 13 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv11')
        trellis = poly2trellis(4,[17 15 13 13 11 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    elseif strcmp(code_mode, 'conv12')
        trellis = poly2trellis(4,[17 15 15 13 13 11 11 7 5 3 1 1]);
        termination = [0, 0, 0];
        tx_code = convenc([msg', termination], trellis);
    end
    tx_code = tx_code';
end