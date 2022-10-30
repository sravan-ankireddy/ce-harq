cd ..; startup; cd conv;

global_settings = 1;

% sim params
nOut = 1;
nMiniFrames = 1000;

nMinFerr = 500;

nFrames = nOut*nMiniFrames;

max_rounds = 10;

% Code parameters
targetCodeRate = 3/4;

N = 400;
K = round(N*targetCodeRate);
R = targetCodeRate;
combining_scheme = "CC";
dec_type = "hard";

if (R == 1/2)
    if (max_rounds == 10)
        SNRdB_low = -10;
        SNRdB_high = -4;
    elseif (max_rounds == 4)
        SNRdB_low = -6;
        SNRdB_high = 0;
    end
elseif (R == 1/3)
    SNRdB_low = -12;
    SNRdB_high = -6;
elseif (R == 1/4)
    SNRdB_low = -16;
    SNRdB_high = -8;
elseif (R == 3/4)
    if (max_rounds == 10)
        SNRdB_low = -8;
        SNRdB_high = -2;
    elseif (max_rounds == 4)
        SNRdB_low = -4;
        SNRdB_high = 2;
    end
elseif (R == 5/6)
    if (max_rounds == 10)
        SNRdB_low = -8;
        SNRdB_high = 2;
    elseif (max_rounds == 4)
        SNRdB_low = -4;
        SNRdB_high = 6;
    end
end
% SNRdB_low = -2;
% SNRdB_high = SNRdB_high - 4;

SNRdB_step = 0.2;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;

channel = "awgn";

if (dec_type == "hard")
    SNRdB_vec = SNRdB_vec + 3;
end

if (channel == "rayleigh")
    SNRdB_vec = SNRdB_vec + 9;
end

num_SNRdB = length(SNRdB_vec);

modulation = 'BPSK';
mod_approx = 0;
comm_mod = 1;

err_thr_ada_scheme = "opt";
channel = "awgn";

res_folder_prefix = 'bler_data';

res_folder_fb = [res_folder_prefix sprintf('/%s/%d/fb/%s/%d',channel, N, modulation,nFrames)];

if ~exist(res_folder_fb,'dir')
    mkdir(res_folder_fb);
end

res_folder_all = res_folder_fb + "/all";
if ~exist(res_folder_all,'dir')
    mkdir(res_folder_all);
end

process_data_fb = 1;

% generate estimates for err_thr
err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));

% run plain harq
err_thr = 0.00;
err_thr_ada_list = err_thr*ones(size(SNRdB_vec));

run_fb;

BLER_HARQ = BLER_vec_FB;
Avg_rounds_HARQ = Avg_rounds_FB;


% run fb with est err_thr
err_thr_ada_scheme = "est";
% if (err_thr_ada_scheme == "est")
%     res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%s/%d/%s/harq_vs_fb_est/%s/%d',channel, N, dec_type, modulation, nFrames)];
% else
%     res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%s/%d/%s/harq_vs_fb_opt/%s/%d',channel, N, dec_type, modulation, nFrames)];
% end

res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%s/%d/%s/harq_vs_fb_est_opt/%s/%d',channel, N, dec_type, modulation, nFrames)];

if ~exist(res_folder_harq_vs_fb,'dir')
    mkdir(res_folder_harq_vs_fb);
end

if (err_thr_ada_scheme == "est")
    min_bler = 1e-5;
    acomp_table = load('lut_data/acomp_400_ns_100000.mat');
    harq_data_path = sprintf('bler_data/awgn/400/hard/harq/BPSK/10000/harq_data_Conv_%d_rate_0.833_rate_0.083_max_rounds_%d.mat',N,max_rounds);
    data = load(harq_data_path);
    for i_rr = 1:size(err_thr_ada_list_est,1)
        remRounds = i_rr;
        for i_ada = 1:size(SNRdB_vec,2)
            err_thr_ada_list_est(i_rr,i_ada) =  err_thr_select(data,acomp_table,targetCodeRate,SNRdB_vec(i_ada),remRounds,min_bler);
        end
    end
end

% err_thr_ada_list_est(1,:) = 0.05;
% err_thr_ada_list_est(2,:) = 0.05;
% err_thr_ada_list_est(3,:) = 0.075;
% err_thr_ada_list_est(4,:) = 0.075;

run_fb;

BLER_FB_est = BLER_vec_FB;
Avg_rounds_FB_est = Avg_rounds_FB;

% run with err thr from grid search

err_thr_grid = 0.00:0.005:0.1;
gs_size = length(err_thr_grid);
err_thr_ada_scheme = "opt";

nFrames_ref = 10000;
res_folder_fb = [res_folder_prefix sprintf('/%s/%d/%s/fb/%s/%d',channel, N, dec_type, modulation, nFrames_ref)];
data_file_name_gs = [res_folder_fb sprintf('/fb_data_Conv_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d.mat', N,R, err_thr_grid(1),err_thr_grid(end), max_rounds)];
gs_data = load(data_file_name_gs);

opt_thr = process_bler_data(gs_data);

err_thr_ada_list = opt_thr.err_thr_opt;
err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));
run_fb;

BLER_FB_opt = BLER_vec_FB;
Avg_rounds_FB_opt = Avg_rounds_FB;

%% Plots

% set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

f = semilogy(SNRdB_vec,BLER_HARQ,'b-o');

hold on;
semilogy(SNRdB_vec,BLER_FB_est,'g-s');
semilogy(SNRdB_vec,BLER_FB_opt,'r-d');

fs = 12;
xlabel('SNR','FontSize',fs);
ylabel('BER','FontSize',fs);

codeRate = R;
leg_HARQ = sprintf('HARQ-%s BLER Rate %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
leg_FB_est = sprintf('FB-%s BLER Rate %.3f, max. %d rounds err thr : est',combining_scheme, codeRate, max_rounds);
leg_FB_opt = sprintf('FB-%s BLER Rate %.3f, max. %d rounds err thr : opt',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB_est, leg_FB_opt, 'Location','southwest','FontSize',fs);

% title_name = sprintf('FB-%s : BLER Conv %d, mod. %s, Rate %.3f, errthr : grid search, max. rounds %d',combining_scheme, N, modulation, R, max_rounds);
title_name = sprintf('FB-%s : BLER Conv %d, mod. %s, Rate %.3f, errthr : grid search and est, max. rounds %d',combining_scheme, N, modulation, R, max_rounds);
title(title_name,'FontSize',fs);

% common_str = [res_folder_harq_vs_fb sprintf('/BLER_HARQ_vs_FB_Conv_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d', N, R, err_thr_ada_scheme, max_rounds)];
common_str = [res_folder_harq_vs_fb sprintf('/BLER_HARQ_vs_FB_Conv_%d_rate_%.3f_err_thr_gs_est_opt_max_rounds_%d', N, R, max_rounds)];
fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);

% set the figure properties
figure('Renderer','painters','Position',[1000 400 800 500]);

f = semilogy(SNRdB_vec,Avg_rounds_HARQ,'b-o');

hold on;
semilogy(SNRdB_vec,Avg_rounds_FB_est,'g-s');
semilogy(SNRdB_vec,Avg_rounds_FB_opt,'r-d');

fs = 12;
xlabel('SNR','FontSize',fs);
ylabel('BER','FontSize',fs);

codeRate = R;
leg_HARQ = sprintf('HARQ-%s Avg Rounds %.3f, max. %d rounds', combining_scheme, codeRate, max_rounds);
leg_FB = sprintf('FB-%s Avg Rounds %.3f, max. %d rounds',combining_scheme, codeRate, max_rounds);

legend(leg_HARQ, leg_FB_est, leg_FB_opt, 'Location','southwest','FontSize',fs);

% title_name = sprintf('FB-%s : Avg rounds Conv %d, mod. %s, Rate %.3f, errthr : grid search %s , max. rounds %d',combining_scheme, N, modulation, R, err_thr_ada_scheme, max_rounds);
title_name = sprintf('FB-%s : Avg rounds Conv %d, mod. %s, Rate %.3f, errthr : grid search and est , max. rounds %d',combining_scheme, N, modulation, R, max_rounds);
title(title_name,'FontSize',fs);

% common_str = [res_folder_harq_vs_fb sprintf('/Avg_rounds_HARQ_vs_FB_Conv_%d_rate_%.3f_err_thr_gs_%s_max_rounds_%d', N, R, err_thr_ada_scheme, max_rounds)];
common_str = [res_folder_harq_vs_fb sprintf('/Avg_rounds_HARQ_vs_FB_Conv_%d_rate_%.3f_err_thr_gs_est_opt_max_rounds_%d', N, R, max_rounds)];

fig_name = common_str + ".fig";
savefig(fig_name);
png_name = common_str + ".png";
saveas(f,png_name);