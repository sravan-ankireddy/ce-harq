global_settings = 1;

% run plain harq
err_thr = 0.00;
err_thr_ada_list = err_thr*ones(size(SNRdB_vec));
err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));

% call this script for running all simulations 
run_fb;

err_thr_ada_scheme = "opt";

BLER_HARQ = BLER_vec_FB;
BER_HARQ = BER_vec_FB;
Avg_rounds_HARQ = Avg_rounds_FB;

% run fb with opt err_thr
err_thr_ada_list = opt_thr.err_thr_opt;
if (inf_rounds == 1)
    err_thr_ada_list = ar_thr.err_thr_opt;
end
run_fb;

BLER_FB_opt = BLER_vec_FB;
BER_FB_opt = BER_vec_FB;
Avg_rounds_FB_opt = Avg_rounds_FB;

% set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

f = semilogy(SNRdB_vec,BLER_HARQ,'b-o');

hold on;
semilogy(SNRdB_vec,BLER_FB_opt,'r-d');

xlabel('SNR');
ylabel('BLER');
grid on;

codeRate = R;
leg_HARQ = sprintf('HARQ');
leg_FB = sprintf('CE-HARQ');
% leg_HARQ = sprintf('HARQ-%s BLER Rate %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
% leg_FB = sprintf('FB-%s BLER Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB);

title_name = sprintf('FB-%s : BLER %s, len %d, %s, Rate %.3f, errthr : %s, max. rounds %d',code_comb_title, combining_scheme, N, modulation, R,err_thr_ada_scheme,  max_rounds);
title(title_name);

if (gs_size > 1)
    common_str = [res_folder_harq_vs_fb sprintf('/BLER_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d',code_comb_str, N, R, err_thr_ada_scheme, max_rounds)];
else
    common_str = [res_folder_harq_vs_fb sprintf('/BLER_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_fixed_%.3f_max_rounds_%d',code_comb_str, N, R,  err_thr_grid(1), max_rounds)];
end

fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);


% set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

f = semilogy(SNRdB_vec,BER_HARQ,'b-o');

hold on;
semilogy(SNRdB_vec,BER_FB_opt,'r-d');

xlabel('SNR');
ylabel('BER');
grid on;

codeRate = R;
leg_HARQ = sprintf('HARQ');
leg_FB = sprintf('CE-HARQ');
% leg_HARQ = sprintf('HARQ-%s BER Rate %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
% leg_FB = sprintf('FB-%s BER Rate %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB);

title_name = sprintf('FB-%s : BER %s, len %d, %s, Rate %.3f, errthr : %s, max. rounds %d',code_comb_title, combining_scheme, N, modulation, R,err_thr_ada_scheme,  max_rounds);
title(title_name);

if (gs_size > 1)
    common_str = [res_folder_harq_vs_fb sprintf('/BER_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d',code_comb_str, N, R, err_thr_ada_scheme, max_rounds)];
else
    common_str = [res_folder_harq_vs_fb sprintf('/BER_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_fixed_%.3f_max_rounds_%d',code_comb_str, N, R,  err_thr_grid(1), max_rounds)];
end


fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);

% AR : set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

f = plot(SNRdB_vec,Avg_rounds_HARQ,'b-o');

hold on;
plot(SNRdB_vec,Avg_rounds_FB_opt,'r-d');

xlabel('SNR');
ylabel('Avg. Rounds');
grid on;

% leg_HARQ = sprintf('HARQ-%s Avg Rounds %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
% leg_FB = sprintf('FB-%s Avg Rounds %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB);

title_name = sprintf('FB-%s : Avg rounds %s, len %d, %s, Rate %.3f, errthr : %s , max. rounds %d',combining_scheme, code_comb_title, N, modulation, R, err_thr_ada_scheme, max_rounds);
title(title_name);

if (gs_size > 1)
    common_str = [res_folder_harq_vs_fb sprintf('/Avg_rounds_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d',code_comb_str, N, R, err_thr_ada_scheme, max_rounds)];
else
    common_str = [res_folder_harq_vs_fb sprintf('/Avg_rounds_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_fixed_%.3f_max_rounds_%d',code_comb_str, N, R,  err_thr_grid(1), max_rounds)];
end

fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);

% SE : set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

SE_HARQ = (K*k/N) * (1 - BLER_HARQ)./Avg_rounds_HARQ;
SE_FB = (K*k/N) * (1 - BLER_FB_opt)./Avg_rounds_FB;

f = plot(SNRdB_vec,SE_HARQ,'b-o');
hold on;
plot(SNRdB_vec,SE_FB,'r-d');

xlabel('SNR');
ylabel('Spectral Efficiency [Bits/s/Hz]');
grid on;

% leg_HARQ = sprintf('HARQ-%s SE %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
% leg_FB = sprintf('FB-%s SE %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB);

title_name = sprintf('FB-%s : Spectral Efficiency %s, len %d, %s, Rate %.3f, errthr : %s , max. rounds %d',code_comb_title, combining_scheme, N, modulation, R, err_thr_ada_scheme, max_rounds);
title(title_name);

if (gs_size > 1)
    common_str = [res_folder_harq_vs_fb sprintf('/Spectral_Efficiency_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d',code_comb_str, N, R, err_thr_ada_scheme, max_rounds)];
else
    common_str = [res_folder_harq_vs_fb sprintf('/Spectral_Efficiency_HARQ_vs_FB_%s_%d_rate_%.3f_err_thr_fixed_%.3f_max_rounds_%d',code_comb_str, N, R,  err_thr_grid(1), max_rounds)];
end

fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);