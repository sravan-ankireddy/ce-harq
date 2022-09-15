function [output, rr_len] = conv_enc(data_in, rate)
    
%     enc = lteConvolutionalEncode(data_in);
%     outlen = ceil(length(data_in)/rate);
%     
%     rr_len = length(enc);
%     output = lteRateMatchConvolutional(enc,outlen);
%     output = double(output);
%       trellis = poly2trellis([5 4],[23 35 0; 0 5 13]);
%       termination = [0, 0, 0, 0, 0, 0, 0, 0];
%       output = convenc([data_in; termination'], trellis);
%       rr_len = length(output);
    trellis = poly2trellis(4,[17 13]);
    termination = [0, 0, 0];
    output = convenc([data_in; termination'], trellis);
    rr_len = length(output);
end