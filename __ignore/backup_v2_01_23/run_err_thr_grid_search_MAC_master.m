startup;

global_settings = 1;
run_grid_search = 1;

err_thr_grid = 0.00:0.01:0.1;
gs_size = length(err_thr_grid);

% sim params
nOut = 1;
nMiniFrames = 500;

nMinFerr = 500;

nFrames = nOut*nMiniFrames;

max_rounds = 4;

% Code parameters
targetCodeRate = 1/2;

N = 800;
PHY_code = "Conv"; % no_code/Conv/LDPC
MAC_code = "Conv"; % no_code/Conv/LDPC

if (MAC_code == "no_code")
    run_grid_search = 0;
end

feedback_mode = "only_PHY"; % MAC_PHY/only_PHY
K = round(N*targetCodeRate);
R = targetCodeRate;
combining_scheme = "CC";
dec_type = "hard";

modulation = 'BPSK';
k = bits_per_symbol(modulation);
M = 2^k;

% FIX ME
ncb = 1;
Nref = 25344;
rvSeq = zeros(1,max_rounds);
max_iter = 6;
nlayers = 1;

% LDPC settings
if (PHY_code == "LDPC")
    dec_type = "unquant";
    targetCodeRate = 1/2;
    K = round(N*targetCodeRate);
    R = targetCodeRate;
    err_thr_grid = 0.02:0.01:0.2;
    gs_size = length(err_thr_grid);

    modulation = 'BPSK';
    k = bits_per_symbol(modulation);
    M = 2^k;

    % PRB settings
    tarCodeLen = 960;
    NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice
    nPRB = round(tarCodeLen/(k*NREPerPRB)); % Vary this to change the code length

    % use nrTBS to get K,R
    if (modulation == 'BPSK')
        tbs = nrTBS('pi/2-BPSK',nlayers,nPRB,NREPerPRB,targetCodeRate); %no. bits in transportBlock
    else
        tbs = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetCodeRate); %no. bits in transportBlock
    end

    % CodeLen based on PRB settings and modulation
    N = nPRB*NREPerPRB*k;
    K = tbs;
    R = K/N;

    % Choose the combining scheme 
    combining_scheme = "CC";
    if (combining_scheme == "IR")
        rvSeq = [0 2 3 1 0 2 3 1 0 2 3 1 0 2 3 1 0 2 3 1 0 2 3 1];
        rvSeq = rvSeq(1:max_rounds);
    else
        rvSeq = zeros(1,max_rounds);
    end
    % SNRdB_low = -4;
    % SNRdB_high = 2;
    if (targetCodeRate == 1/3)
        if (max_rounds == 10)
            SNRdB_low = -10;
            SNRdB_high = -4;
        elseif (max_rounds == 4)
            SNRdB_low = -8;
            SNRdB_high = -4;
        end
    elseif (targetCodeRate == 1/2)
        if (max_rounds == 4)
            SNRdB_low = -8;
            SNRdB_high = -2;
        end   
    elseif (targetCodeRate == 3/4)
        if (max_rounds == 4)
            SNRdB_low = -6;
            SNRdB_high = 0;
        end
    elseif (targetCodeRate == 0.9)
        if (max_rounds == 4)
            SNRdB_low = -4;
            SNRdB_high = 2;
        end
    end
    rv = 0;
elseif (PHY_code == "Conv")
    if (R == 1/2)
        if (max_rounds == 10)
            SNRdB_low = -10;
            SNRdB_high = -4;
        elseif (max_rounds == 4)
            SNRdB_low = -6;
            SNRdB_high = 0;
        end
    elseif (R == 1/3)
        SNRdB_low = -9;
        SNRdB_high = -3;
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
elseif (PHY_code == "no_code")
    targetCodeRate = 1;
    K = N;
    R = 1;
    if (max_rounds == 4)
        SNRdB_low = 0;
        SNRdB_high = 10;
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

if (modulation == "QPSK")
    SNRdB_vec = SNRdB_vec + 3;
end
num_SNRdB = length(SNRdB_vec);

% Error Stats
BER_vec_FB_gs = zeros(gs_size,num_SNRdB);
BLER_vec_FB_gs = zeros(gs_size,num_SNRdB);

BER_vec_pr_FB_gs = zeros(gs_size,max_rounds,num_SNRdB);
BLER_vec_pr_FB_gs = zeros(gs_size,max_rounds,num_SNRdB);

Avg_rounds_FB_gs = zeros(gs_size,num_SNRdB);


mod_approx = 0;
comm_mod = 1;

err_thr_ada_scheme = "opt";

res_folder_prefix = 'bler_data_master';

code_comb = sprintf('/PHY_%s_MAC_%s', PHY_code, MAC_code);

res_folder_fb = [res_folder_prefix sprintf('/%s/%s/%s/%d/%s/fb/%s/%d',feedback_mode, code_comb, channel, N, dec_type, modulation, nFrames)];

if (err_thr_ada_scheme == "est")
    res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%s/%s/%s/%d/%s/harq_vs_fb_est/%s/%d',feedback_mode, code_comb, channel, N, dec_type, modulation,nFrames)];
else
    res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%s/%s/%s/%d/%s/harq_vs_fb_opt/%s/%d',feedback_mode, code_comb, channel, N, dec_type, modulation,nFrames)];
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
        run_fb_MAC_master;

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
    grid on;
    legend(leg_str);
    title_str = sprintf('FB-(%s PHY, %s MAC) scheme : BLER %s %d mod. %s Rate %.3f max.rounds %d',PHY_code, MAC_code, combining_scheme, N, modulation, R, max_rounds);
    title(title_str);
    bler_common_str = [res_folder_fb sprintf('/FB_BLER_%s_%s_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d_numF_%d',PHY_code, MAC_code, N, R, err_thr_grid(1), err_thr_grid(end), max_rounds, nFrames)];
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
    grid on;
    legend(leg_str);
    title_str = sprintf('FB-(%s PHY, %s MAC) scheme : Avg rounds %s %d mod. %s Rate %.3f max. rounds %d',PHY_code, MAC_code, combining_scheme, N, modulation, R, max_rounds);
    title(title_str);
    ar_common_str = [res_folder_fb sprintf('/Avg_rounds_%s_%s_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d_numF_%d',PHY_code, MAC_code, N, R, err_thr_grid(1), err_thr_grid(end), max_rounds, nFrames)];
    filename_AR_fig = ar_common_str + ".fig";
    filename_AR_png = ar_common_str + ".png";

    savefig(filename_AR_fig);
    saveas(f,filename_AR_png);

    % set the figure properties for SE plots
    f = figure('Renderer','painters','Position',[1000 400 800 500]);
    
    leg_str = {};
    for i_e = 1:length(err_thr_grid)
        SE = (K*k/N) *(1 - squeeze(BLER_vec_FB_gs(i_e,:))) ./squeeze(Avg_rounds_FB_gs(i_e,:));
        semilogy(SNRdB_vec,SE);
        hold on;
        leg_str{end+1} = sprintf('Err thr %.3f',err_thr_grid(i_e));
    end
    xlabel('SNR');
    ylabel('Spectral Efficiency [bits/s/Hz]');
    grid on;
    legend(leg_str);
    title_str = sprintf('FB-(%s PHY, %s MAC) scheme : Avg rounds %s %d mod. %s Rate %.3f max. rounds %d',PHY_code, MAC_code, combining_scheme, N, modulation, R, max_rounds);
    title(title_str);
    se_common_str = [res_folder_fb sprintf('/SE_%s_%s_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d_numF_%d',PHY_code, MAC_code, N, R, err_thr_grid(1), err_thr_grid(end), max_rounds, nFrames)];
    filename_SE_fig = se_common_str + ".fig";
    filename_SE_png = se_common_str + ".png";

    savefig(filename_SE_fig);
    saveas(f,filename_SE_png);

    ber_data = BER_vec_pr_FB_gs;
    bler_data = BLER_vec_pr_FB_gs;
    ar_data = Avg_rounds_FB_gs;
    snr_data = SNRdB_vec;
end

% Store data
if (run_grid_search == 1)
    outer_str = sprintf('/fb_data_%s_%s_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d.mat', PHY_code, MAC_code, N, R, err_thr_grid(1),err_thr_grid(end), max_rounds);
    data_file_name_gs = [res_folder_fb outer_str];
    save(data_file_name_gs,'ber_data','bler_data','ar_data','snr_data','err_thr_grid');
    gs_data = load(data_file_name_gs);
    opt_thr = process_bler_data(gs_data);
else
    if (MAC_code == "no_code")
        opt_thr.err_thr_opt = zeros(1,length(SNRdB_vec));
    else
        gs_data = load(data_file_name_gs);
        opt_thr = process_bler_data(gs_data);
        nFrames_ref = 100000;
        res_folder_fb = [res_folder_prefix sprintf('/%s/%d/%s/fb/%s/%d',channel, N, dec_type, modulation, nFrames_ref)];
        outer_str = sprintf('/fb_data_%s_%s_%d_rate_%.3f_err_thr_%.3f_to_%.3f_max_rounds_%d.mat', PHY_code, MAC_code, N, R, err_thr_grid(1),err_thr_grid(end), max_rounds);
        data_file_name_gs = [res_folder_fb outer_str];
    end
end

run_harq_vs_fb_MAC_master;