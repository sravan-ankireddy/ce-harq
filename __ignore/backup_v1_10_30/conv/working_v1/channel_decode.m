function msg_hat = channel_decode(rx_code, code_mode, interleaver)

    if nargin > 2
        intrlvrIndices = interleaver;
    else
        intrlvrIndices = 0;
    end

    if strcmp(code_mode, 'uncoded')
        msg_hat = rx_code;
    elseif strcmp(code_mode, 'turbo')
        trellis = poly2trellis(4,[13 15],13);
        numiter = 4;
        turbodec = comm.TurboDecoder(trellis,intrlvrIndices,numiter);
        msg_hat = turbodec(2.0*rx_code'-1.0)';
    elseif strcmp(code_mode, 'conv2')
        trellis = poly2trellis(4,[17 13]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv3')
        trellis = poly2trellis(4,[17 15 13]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv4')
        trellis = poly2trellis(4,[17 15 13 11]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv5')
        trellis = poly2trellis(4,[17 15 13 11 7]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv6')
        trellis = poly2trellis(4,[17 15 13 11 7 5]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv7')
        trellis = poly2trellis(4,[17 15 13 11 7 5 3]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv8')
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv9')
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded; 
    elseif strcmp(code_mode, 'conv10')
        trellis = poly2trellis(4,[17 15 13 13 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv11')
        trellis = poly2trellis(4,[17 15 13 13 11 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif strcmp(code_mode, 'conv12')
        trellis = poly2trellis(4,[17 15 15 13 13 11 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        msg_hat = decoded;  
    end
end