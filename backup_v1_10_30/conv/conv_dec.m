function [msg_hat] = conv_dec(llr, rate, dec_type)

    % quatize the llr to get hard est
    msg_hat = 0;
    rx_code = llr;

    if (nargin < 3)
        dec_type = 'hard';
    end
   
    if (rate == 5/6)
        trellis = poly2trellis(4,[17 13]);
        tb_len = min(15, round(length(llr)/3));
        puncpat = [1;1;0;1;0];
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type, puncpat);
        % FIX ME
        if (length(rx_code) <= 408)
            decoded = decoded(1:end-7);
        else
            decoded = decoded(1:end-4);
        end
        msg_hat = decoded;
    elseif (rate == 3/4)
        trellis = poly2trellis(4,[17 13]);
        tb_len = min(15, round(length(llr)/3));
        puncpat = [1;1;0];
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type, puncpat);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/2)
        trellis = poly2trellis(4,[17 13]);
        tb_len = min(15, round(length(llr)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/3)
        trellis = poly2trellis(4,[17 15 13]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/4)
        trellis = poly2trellis(4,[17 15 13 11]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/5)
        trellis = poly2trellis(4,[17 15 13 11 7]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/6)
        trellis = poly2trellis(4,[17 15 13 11 7 5]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/7)
        trellis = poly2trellis(4,[17 15 13 11 7 5 3]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/8)
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/9)
        trellis = poly2trellis(4,[17 15 13 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/10)
        trellis = poly2trellis(4,[17 15 13 13 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/11)
        trellis = poly2trellis(4,[17 15 13 13 11 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    elseif (rate == 1/12)
        trellis = poly2trellis(4,[17 15 15 13 13 11 11 7 5 3 1 1]);
        tb_len = min(15, round(length(rx_code)/3));
        decoded = vitdec(rx_code, trellis, tb_len, 'term', dec_type);
        decoded = decoded(1:end-3);
        msg_hat = decoded;
    end
end