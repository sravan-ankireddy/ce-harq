
% BSC (p)
p_vec = 0:0.001:0.01;
x_err_vec = zeros(size(p_vec));
BER_vec = zeros(size(p_vec));
n_samples = 1e4;
spa = 0.05;

for i_p = 1:length(p_vec)
    p = p_vec(i_p);

    % Create Huffman dictionary
    H = 6; % Huffman dictionary size

    symbols = 1:6; 
    p_sym = [.5 .125 .125 .125 .0625 .0625];
    
    dict = huffmandict(symbols,p_sym);
    
    inputSig = randsrc(100,1,[symbols;p_sym]);
    
    code = huffmanenco(inputSig,dict);

    % pass through bsc(p)
    y = bsc(code,p);
    
    sig = huffmandeco(y,dict);

    binarySig = de2bi(inputSig);
    binaryOutputSig = de2bi(sig);

    BER_vec(i_p) = mean(binarySig~=binaryOutputSig,"all");

    disp(BER_vec(i_p));

%     isequal(inputSig,sig)


%     % Create a sparse vector and Compress using Huffman coding
%     sym = [00 01 10 11];
%     p_sym = [(1-spa)^2 spa*(1-spa) spa*(1-spa) spa^2];
% 
%     dict = huffmandict(sym,p_sym);
% 
%     u = randsrc(100,1,[sym;p_sym]);
% 
%     x = huffmanenco(u,dict);
% 
%     % pass through bsc(p)
%     y = bsc(x,p);
% 
%     % Decode using Huffman
%     x_hat = huffmandeco(y,dict);
end

plot(p_vec,BER_vec);