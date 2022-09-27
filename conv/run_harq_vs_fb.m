global_settings = 1;

gs_data = load(data_file_name_gs);

opt_thr = process_bler_data(gs_data);

% run plain harq
err_thr = 0.00;
err_thr_ada_list = err_thr*ones(size(SNRdB_vec));
err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));

err_thr_ada_list_est(1,:) = 0.05;
err_thr_ada_list_est(2,:) = 0.06;
err_thr_ada_list_est(3,:) = 0.075;
err_thr_ada_list_est(4,:) = 0.075;
err_thr_ada_scheme = "est";
run_fb;

err_thr_ada_scheme = "opt";

BLER_HARQ = BLER_vec_FB;
Avg_rounds_HARQ = Avg_rounds_FB;

% run fb with opt err_thr
err_thr_ada_list = opt_thr.err_thr_opt;
run_fb;

BLER_FB_opt = BLER_vec_FB;
Avg_rounds_FB_opt = Avg_rounds_FB;

% set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

f = semilogy(SNRdB_vec,BLER_HARQ,'b-o');

hold on;
semilogy(SNRdB_vec,BLER_FB_opt,'r-d');

fs = 12;
xlabel('SNR','FontSize',fs);
ylabel('BLER','FontSize',fs);

codeRate = R;
leg_HARQ = sprintf('FB-HARQ BLER Rate %.3f, max. %d rounds', codeRate, max_rounds);
leg_FB = sprintf('FB-%s BLER Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB, 'Location','southwest','FontSize',fs);

title_name = sprintf('FB-%s : BLER Conv %d, mod. %s, Rate %.3f, errthr : grid search %s, max. rounds %d',combining_scheme, N, modulation, R,err_thr_ada_scheme,  max_rounds);
title(title_name,'FontSize',fs);

common_str = [res_folder_harq_vs_fb sprintf('/BLER_HARQ_vs_FB_Conv_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d', N, R, err_thr_ada_scheme, max_rounds)];

fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);

% set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

f = semilogy(SNRdB_vec,Avg_rounds_HARQ,'b-o');

hold on;
semilogy(SNRdB_vec,Avg_rounds_FB_opt,'r-d');

fs = 12;
xlabel('SNR','FontSize',fs);
ylabel('Avg. Rounds','FontSize',fs);

codeRate = R;
leg_HARQ = sprintf('FB-HARQ Avg Rounds %.3f, max. %d rounds', codeRate, max_rounds);
leg_FB = sprintf('FB-%s Avg Rounds %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB, 'Location','southwest','FontSize',fs);

title_name = sprintf('FB-%s : Avg rounds Conv %d, mod. %s, Rate %.3f, errthr : grid search %s , max. rounds %d',combining_scheme, N, modulation, R, err_thr_ada_scheme, max_rounds);
title(title_name,'FontSize',fs);

common_str = [res_folder_harq_vs_fb sprintf('/Avg_rounds_HARQ_vs_FB_Conv_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d', N,R, err_thr_ada_scheme, max_rounds)];

fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);