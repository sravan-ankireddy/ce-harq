function [output] = conv_dec(llr, block_len)

%     raterec = lteRateRecoverConvolutional(llr,block_len);
% 
%     decBits = lteConvolutionalDecode(raterec);

%     output = double(decBits);
    llr = llr > 0;
%     trellis = poly2trellis([5 4],[23 35 0; 0 5 13]);
%     tb_len = min(15, round(length(llr)/3));
%     decoded = vitdec(llr, trellis, tb_len, 'term', 'hard');
%     decoded = decoded(1:end-3);
%     output = decoded;

    trellis = poly2trellis(4,[17 13]);
        tb_len = min(15, round(length(llr)/3));
        decoded = vitdec(llr, trellis, tb_len, 'term', 'hard');
        decoded = decoded(1:end-3);
        output = decoded;
end