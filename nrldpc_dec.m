function [output, crc] = nrldpc_dec(chOut, rate, modulation, block_len, max_iter, rv, nlayers, bgn)
    CRC_L = 16;

    raterec = nrRateRecoverLDPC(chOut, block_len, rate, rv, modulation, nlayers);

    decBits = nrLDPCDecode(raterec, bgn, max_iter, 'Algorithm','Normalized min-sum');

    % Code block desegmentation and CRC decoding
    [blk, ~] = nrCodeBlockDesegmentLDPC(decBits, bgn, block_len+CRC_L);

    % Transport block CRC decoding
    [output, crc] = nrCRCDecode(blk,  '16');
    output = double(output);
end