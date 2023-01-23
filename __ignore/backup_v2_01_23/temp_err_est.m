% max_rounds = 4;
% SNRdB_vec = 1:0.2:9;
err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));
targetCodeRate = 3/4;
% if (err_thr_ada_scheme == "est")
    min_bler = 1e-5;
    acomp_table = load('lut_data/acomp_400_ns_100000.mat');
    data = load(harq_data_path);
    for i_rr = 1:size(err_thr_ada_list_est,1)
        remRounds = i_rr;
        for i_ada = 1:size(SNRdB_vec,2)
            err_thr_ada_list_est(i_rr,i_ada) =  err_thr_select(data,acomp_table,targetCodeRate,SNRdB_vec(i_ada),remRounds,min_bler);
        end
    end
% end
disp(err_thr_ada_list_est);

opt_thr = process_bler_data(gs_data);

err_thr_ada_list = opt_thr.err_thr_opt;