global_settings = 1;

unit_test = 1;

if ~unit_test
    nPRB_list = [8, 20, 60];
    targetCodeRate_list = [0.8, 0.9, 0.2];
    max_iter = 6;
    err_thr_list = [0.02, 0.05];
    num_blocks = 10e0;
else
    nPRB_list = [8];
    targetCodeRate_list = [0.8];
    max_iter = 6;
    err_thr_list = [0.02];
    num_blocks = 10e2;
end

for i_n = 1:length(nPRB_list)
    for i_t = 1:length(targetCodeRate_list)
        for i_e = 1:length(err_thr_list)
            nPRB = nPRB_list(i_n);
            targetCodeRate = targetCodeRate_list(i_t);
            err_thr = err_thr_list(i_e);
            % Call the main script
            harq_vs_fb;
        end
    end
end

clear;