cd ..; startup; cd conv;

global_settings = 1;
run_grid_search = 1;

err_thr_grid = 0.00:0.01:0.15;
gs_size = length(err_thr_grid);

% sim params
nOut = 10;
nMiniFrames = 50;

nFrames = nOut*nMiniFrames;

max_rounds = 4;

% Code parameters
targetCodeRate = 3/4;

N = 400;
K = round(N*targetCodeRate);
R = targetCodeRate;
combining_scheme = "CC";
dec_type = "unquant";

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
end

SNRdB_low = SNRdB_low - 4;
SNRdB_high = SNRdB_high - 4;

SNRdB_step = 0.2;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;

num_SNRdB = length(SNRdB_vec);

% Error Stats
BER_vec_FB_gs = zeros(gs_size,num_SNRdB);
BLER_vec_FB_gs = zeros(gs_size,num_SNRdB);

BER_vec_pr_FB_gs = zeros(gs_size,max_rounds,num_SNRdB);
BLER_vec_pr_FB_gs = zeros(gs_size,max_rounds,num_SNRdB);

Avg_rounds_FB_gs = zeros(gs_size,num_SNRdB);


modulation = 'BPSK';
mod_approx = 0;
comm_mod = 1;

err_thr_ada_scheme = "opt";
channel = "awgn";

res_folder_prefix = 'bler_data';

res_folder_fb = [res_folder_prefix sprintf('/%d/%s/fb/%s/%d', N, dec_type, modulation, nFrames)];

if (err_thr_ada_scheme == "est")
    res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%d/%s/harq_vs_fb_est/%s/%d',N, dec_type, modulation,nFrames)];
else
    res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%d/%s/harq_vs_fb_opt/%s/%d',N, dec_type, modulation,nFrames)];
end

if ~exist(res_folder_fb,'dir')
    mkdir(res_folder_fb);
end

if ~exist(res_folder_harq_vs_fb,'dir')
    mkdir(res_folder_harq_vs_fb);
end

res_folder_all = res_folder_fb + "/all";
if ~exist(res_folder_all,'dir')
    mkdir(res_folder_all);
end

process_data_fb = 1;

if (run_grid_search == 1)
    for i_e = 1:gs_size
        err_thr = err_thr_grid(i_e);
        err_thr_ada_list = err_thr*ones(size(SNRdB_vec));
        err_thr_ada_list_est = 0;
        run_fb;

        % Error Stats
        BER_vec_FB_gs(i_e,:) = BER_vec_FB;
        BLER_vec_FB_gs(i_e,:) = BLER_vec_FB;
        
        BER_vec_pr_FB_gs(i_e,:,:) = BER_vec_pr_FB;
        BLER_vec_pr_FB_gs(i_e,:,:) = BLER_vec_pr_FB;
        
        Avg_rounds_FB_gs(i_e,:) = Avg_rounds_FB;
    end

    %% Plot the cumulative results
    % set the figure properties for BLER plots
    f = figure('Renderer','painters','Position',[1000 400 800 500]);

    leg_str = {};
    for i_e = 1:length(err_thr_grid)
        semilogy(SNRdB_vec,squeeze(BLER_vec_FB_gs(i_e,:)));
        hold on;
        leg_str{end+1} = sprintf('Err thr %.3f',err_thr_grid(i_e));
    end
    xlabel('SNR');
    ylabel('BLER');
    legend(leg_str);
    title_str = sprintf('FB-%s scheme : BLER LDPC %d mod. %s Rate %.3f max.rounds %d',combining_scheme, N, modulation, R, max_rounds);
    title(title_str);
    bler_common_str = [res_folder_fb sprintf('/FB_BLER_LDPC_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d_numF_%d',N,R, err_thr_grid(1), err_thr_grid(end), max_rounds, nFrames)];
    filename_BLER_fig = bler_common_str + ".fig";
    filename_BLER_png = bler_common_str + ".png";

    savefig(filename_BLER_fig);
    saveas(f,filename_BLER_png);

    % set the figure properties for AR plots
    f = figure('Renderer','painters','Position',[1000 400 800 500]);

    leg_str = {};
    for i_e = 1:length(err_thr_grid)
        semilogy(SNRdB_vec,squeeze(Avg_rounds_FB_gs(i_e,:)));
        hold on;
        leg_str{end+1} = sprintf('Err thr %.3f',err_thr_grid(i_e));
    end
    xlabel('SNR');
    ylabel('Avg. rounds');
    legend(leg_str);
    title_str = sprintf('FB-%s scheme : AR LDPC %d mod. %s Rate %.3f max. rounds %d',combining_scheme, N, modulation, R, max_rounds);
    title(title_str);
    ar_common_str = [res_folder_fb sprintf('/FB_AR_LDPC_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d_numF_%d',N,R, err_thr_grid(1), err_thr_grid(end), max_rounds, nFrames)];
    filename_AR_fig = ar_common_str + ".fig";
    filename_AR_png = ar_common_str + ".png";

    savefig(filename_AR_fig);
    saveas(f,filename_AR_png);

    ber_data = BER_vec_pr_FB_gs;
    bler_data = BLER_vec_pr_FB_gs;
    ar_data = Avg_rounds_FB_gs;
    snr_data = SNRdB_vec;
end

% Store data
if (run_grid_search == 1)
    data_file_name_gs = [res_folder_fb sprintf('/fb_data_Conv_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d.mat', N,R, err_thr_grid(1),err_thr_grid(end), max_rounds)];
    save(data_file_name_gs,'ber_data','bler_data','ar_data','snr_data','err_thr_grid');
else
    nFrames_ref = 10000
    res_folder_fb = [res_folder_prefix sprintf('/%d/fb/%s/%d',N, modulation,nFrames_ref)];
    data_file_name_gs = [res_folder_fb sprintf('/fb_data_Conv_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d.mat', N,R, err_thr_grid(1),err_thr_grid(end), max_rounds)];
end

run_harq_vs_fb;