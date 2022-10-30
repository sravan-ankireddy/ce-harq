function [Rc_vec, sp_vec] = generate_acomp_table(N)
    % Input sparsity
    sp_vec = 0.5:0.5:10;
    
    len_vec = zeros(1,length(sp_vec));
    
    num_samples = 10e4;
    tic;
    for i_s = 1:length(sp_vec)
        s = sp_vec(i_s);
        disp(i_s);
        temp_len = 0;
        parfor i_n = 1:num_samples
        
            % Generate random inputs of length N
            x0 = rand(N,1);
            x = x0;
            x(x0>0.01*s) = 0;
            x(x0<=0.01*s) = 1;
        
            % Define the compressor
            counts = [100-round(s) round(s)];
        
            x_comp = arithenco(x+1,counts);
        
            temp_len = temp_len + length(x_comp);
        end

        len_vec(i_s) = len_vec(i_s) + temp_len;
        
        len_vec(i_s) = len_vec(i_s)/num_samples;
    end
    toc;
    Rc_vec = len_vec/(N);
    sp_vec = sp_vec/100;
    % Save to a mat file
    filename = sprintf('lut_data/acomp_%d_ns_%d.mat',N,num_samples);
    save(filename,'Rc_vec','sp_vec');
end