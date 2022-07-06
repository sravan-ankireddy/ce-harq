% Input sparsity
sp_vec = 1:10;

% Length of input
N = 1000;

len_vec = zeros(1,length(sp_vec));

num_samples = 10e2;
tic;
for i_s = 1:length(sp_vec)
    s = sp_vec(i_s);
    disp(i_s);
    for i_n = 1:num_samples
    
        % Generate random inputs of length N
        x0 = rand(N,1);
        x = x0;
        x(x0>0.01*s) = 0;
        x(x0<=0.01*s) = 1;
    
        % Define the compressor
        counts = [100-s s];
    
        x_comp = arithenco(x+1,counts);
    
        len_vec(i_s) = len_vec(i_s) + length(x_comp);
    end
    
    len_vec(i_s) = len_vec(i_s)/num_samples;
end
toc;
Avg_compression = len_vec/(N);