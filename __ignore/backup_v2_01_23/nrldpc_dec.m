function [output, crc] = nrldpc_dec(raterec, block_len, max_iter, bgn)

    CRC_L = 16;
    decBits = nrLDPCDecode(raterec, bgn, max_iter, 'Algorithm','Normalized min-sum');

    % Code block desegmentation and CRC decoding
    [blk, ~] = nrCodeBlockDesegmentLDPC(decBits, bgn, block_len+CRC_L);

    % Transport block CRC decoding
    [output, crc] = nrCRCDecode(blk,  '16');
    output = double(output);
end