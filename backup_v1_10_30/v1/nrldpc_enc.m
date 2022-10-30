function [output] = nrldpc_enc(data_in, rate, modulation, rv, bgn, nlayers)
    in = int8(data_in);
    tbIn = nrCRCEncode(in,  '16');
    cbsIn = nrCodeBlockSegmentLDPC(tbIn, bgn);
    enc = nrLDPCEncode(cbsIn, bgn);
    outlen = ceil(length(data_in)/rate);
    output = nrRateMatchLDPC(enc,outlen,rv,modulation,nlayers);
    output = double(output);
end