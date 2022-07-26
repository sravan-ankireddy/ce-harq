close all; clear all; clc;

% Figure settings
startup;

global_settings = 1;

unit_test = 0;

% run params
run_er_thr_grid_search = 1;

% use either estimate or opt values from LUT
err_thr_ada_scheme = "opt";
targetCodeRate = 0.2;

% generate acomp LUT
acomp_lut_path = "lut_data/acomp_960_ns_100000.mat";
if ~isfile(acomp_lut_path)
   acomp_table = generate_acomp_table;
else
    acomp_table = load(acomp_lut_path);
end

% constant settings
% FIX ME
ncb = 1;
Nref = 25344;
max_iter = 3;
mod_approx = 0;
qam_mod = 0;
max_rounds = 10;
nFrames = 10e2;
min_bler = 1e-4; % URLLC requirement 1e-5
modulation = 'QPSK';

% PRB settings
nPRB = 20; % Vary this to change the code length
nlayers = 1;
NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice

% CodeLen based on PRB settings
N = nPRB*NREPerPRB;

% use nrTBS to get K,R
tbs = nrTBS(modulation,nlayers,nPRB,NREPerPRB,targetCodeRate);

% Use generated code params
M = bits_per_symbol(modulation);
K = tbs/M;
R = K/N;

% Choose the combining scheme 
combining_scheme = "IR";
if (combining_scheme == "IR")
    rvSeq = [0 2 3 1 0 2 3 1 0 2 3 1 0 2 3 1 0 2 3 1 0 2 3 1];
    rvSeq = rvSeq(1:max_rounds);
else
    rvSeq = zeros(1,max_rounds);
end

% path to results folder
if (combining_scheme == "IR")
    res_folder_prefix = 'bler_data_ir';
else
    res_folder_prefix = 'bler_data_cc';
end

res_folder_fb = [res_folder_prefix sprintf('/%d/fb/%d',N,nFrames)];
if (err_thr_ada_scheme == "est")
    res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%d/harq_vs_fb_est/%d',N,nFrames)];
else
    res_folder_harq_vs_fb = [res_folder_prefix sprintf('/%d/harq_vs_fb_opt/%d',N,nFrames)];
end

% FB params
res_folder = res_folder_fb;
err_thr_list = 0.00:0.005:0.05;
if (combining_scheme == "IR")
    SNRdB_low = -6;
    SNRdB_high =-2;
    if (max_rounds == 5 && max_iter == 3)
        SNRdB_low = -2;
        SNRdB_high =2;
    end

    if (max_rounds == 5 && max_iter == 6)
        SNRdB_low = -5;
        SNRdB_high =-1;
    end

    if (max_rounds == 3 && max_iter == 6)
        SNRdB_low = -2;
        SNRdB_high = 2;
    end
else
    SNRdB_low = -6;
    SNRdB_high =-1;
end

SNRdB_step = 0.2;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
num_SNRdB = length(SNRdB_vec);
snr_data = SNRdB_vec;

if (targetCodeRate == 0.2)
    SNRdB_vec = SNRdB_vec - 7;
end

if (unit_test)
    err_thr_list = 0.01:0.01:0.02;
    nFrames = 10e1;
end

err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));

% Initiate parpool
parpool_size = min(64,feature('numcores'));
pp = gcp('nocreate'); % If no pool, do not create new one.
if (parpool_size ~= 8 && parpool_size ~= 10 && isempty(pp))
    parpool(parpool_size);
end

%% Feedback grid_search data
if (run_er_thr_grid_search == 1)
       
    if ~exist(res_folder,'dir')
        mkdir(res_folder);
    end
    
    res_folder_all = res_folder + "/all";
    if ~exist(res_folder_all,'dir')
        mkdir(res_folder_all);
    end

    % Store full error data : 2xerr_thr_listxmax_roundsxnum_SNR
    err_data_all_thr = zeros(2,length(err_thr_list),max_rounds,num_SNRdB);

    % Store Avg. rounds data : err_thr_listxnum_SNR
    ar_data_all_thr = zeros(length(err_thr_list),num_SNRdB);
    for i_e = 1:length(err_thr_list)
        err_thr = err_thr_list(i_e);
        err_thr_ada_list = err_thr*ones(size(SNRdB_vec));

        % Call the main script
        process_data_fb = 1;
        run_fb;

        % Save the results            
        data_file_name = [res_folder sprintf('/fb_lut_data_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d.mat',N,R, max_iter, err_thr_list(1), err_thr_list(end), max_rounds, qam_mod, mod_approx, nFrames)];

        save(data_file_name,'err_data_all_thr','ar_data_all_thr','snr_data','err_thr_list');
    end

    %% Plot the cumulative results
    % set the figure properties for BLER plots
    f = figure('Renderer','painters','Position',[1000 400 800 500]);

    leg_str = {};
    for i_e = 1:length(err_thr_list)
        semilogy(SNRdB_vec,squeeze(err_data_all_thr(2,i_e,end,:)));
        hold on;
        leg_str{end+1} = sprintf('Err thr %.3f',err_thr_list(i_e));
    end
    xlabel('SNR');
    ylabel('BLER');
    legend(leg_str);
    title_str = sprintf('FB scheme : BLER LDPC %d Rate %.3f max. iter %d max.rounds %d',N,R, max_iter, max_rounds);
    title(title_str);
    bler_common_str = [res_folder sprintf('/FB_BLER_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d',N,R, max_iter, err_thr_list(1), err_thr_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
    filename_BLER_fig = bler_common_str + ".fig";
    filename_BLER_png = bler_common_str + ".png";

    savefig(filename_BLER_fig);
    saveas(f,filename_BLER_png);

    % set the figure properties for AR plots
    f = figure('Renderer','painters','Position',[1000 400 800 500]);

    leg_str = {};
    for i_e = 1:length(err_thr_list)
        semilogy(SNRdB_vec,squeeze(ar_data_all_thr(i_e,:)));
        hold on;
        leg_str{end+1} = sprintf('Err thr %.3f',err_thr_list(i_e));
    end
    xlabel('SNR');
    ylabel('Avg. rounds');
    legend(leg_str);
    title_str = sprintf('FB scheme : AR LDPC %d Rate %.3f max. iter %d max. rounds %d',N,R, max_iter, max_rounds);
    title(title_str);
    ar_common_str = [res_folder sprintf('/FB_AR_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d',N,R, max_iter, err_thr_list(1), err_thr_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
    filename_AR_fig = ar_common_str + ".fig";
    filename_AR_png = ar_common_str + ".png";

    savefig(filename_AR_fig);
    saveas(f,filename_AR_png);
end

%% Compare HARQ & FB based on grid_search data

% path to gris_search data
actualRate = R;
gs_data_path = [res_folder_fb sprintf('/fb_lut_data_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_to_%.3f_max_rounds_%d_qm_0_ma_0_numF_%d.mat',N,actualRate,max_iter,err_thr_list(1),err_thr_list(end),max_rounds,nFrames)];

res_folder = res_folder_harq_vs_fb;
if ~exist(res_folder,'dir')
    mkdir(res_folder);
end

% generate estimates for err_thr
err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));
if (err_thr_ada_scheme == "est")
    data = load(harq_lut_path);
    for i_rr = 1:size(err_thr_ada_list_est,1)
        remRounds = i_rr;
        for i_ada = 1:size(SNRdB_vec,2)
            err_thr_ada_list_est(i_rr,i_ada) =  err_thr_select(data,acomp_table,targetCodeRate,SNRdB_vec(i_ada),remRounds,min_bler);
        end
    end
end

% HARQ
% Store full error data : 2xmax_roundsxnum_SNR : first column of grid search data with err_thr = 0
data = load(gs_data_path);
err_data_harq = squeeze(data.err_data_all_thr(:,1,:,:));

ar_data_harq = data.ar_data_all_thr(1,:);

% look up opt err_thr for per SNR
[bler_out, ar_out, snr_out] = process_bler_data(data);
err_thr_ada_list = bler_out.err_thr_opt;

% extrapolate at edges to match size of SNRdB_vec
if (snr_out(1) > SNRdB_low)
    nz = (snr_out(1) - SNRdB_low)/SNRdB_step;
    err_thr_ada_list = [err_thr_ada_list(1)*ones(nz,1); err_thr_ada_list];
end
if (snr_out(end) < SNRdB_high)
    nz = (SNRdB_high - snr_out(end))/SNRdB_step;
    err_thr_ada_list = [err_thr_ada_list; err_thr_ada_list(end)*ones(nz,1)];
end

% FB
% Store full error data : 2xmax_roundsxnum_SNR
err_data_fb = zeros(size(err_data_harq));

% call the main script and run new sim with opt/est err_thr
process_data_fb = 0;
run_fb;
err_data_fb(1,:,:) = BER_vec_pr_FB;
err_data_fb(2,:,:) = BLER_vec_pr_FB;
ar_data_fb = Avg_rounds_FB;

% Save the results            
data_file_name = [res_folder sprintf('/harq_vs_fb_data_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_ada_max_rounds_%d_qm_%d_ma_%d_numF_%d.mat',N,R, max_iter, max_rounds, qam_mod, mod_approx, nFrames)];

save(data_file_name,'err_data_fb','ar_data_fb','snr_data','err_thr_ada_list');

%% Plot and compare results

%% Plot the cumulative results
% set the figure properties for BLER plots
f = figure('Renderer','painters','Position',[1000 400 800 500]);

semilogy(SNRdB_vec,squeeze(err_data_harq(2,end,:)));
hold on;
    semilogy(SNRdB_vec,squeeze(err_data_fb(2,end,:)));
xlabel('SNR');
ylabel('BLER');
legend("HARQ","FB");
title_str = sprintf('HARQ vs FB scheme opt. err thr : BLER LDPC %d Rate %.3f max. iter %d max.rounds %d',N,R, max_iter, max_rounds);
title(title_str);
bler_common_str = [res_folder sprintf('/BLER_HARQ_vs_FB_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_ada_max_rounds_%d_qm_%d_ma_%d_numF_%d',N,R, max_iter, max_rounds, qam_mod, mod_approx, nFrames)];
filename_BLER_fig = bler_common_str + ".fig";
filename_BLER_png = bler_common_str + ".png";

savefig(filename_BLER_fig);
saveas(f,filename_BLER_png);

% set the figure properties for BLER plots
f = figure('Renderer','painters','Position',[1000 400 800 500]);

semilogy(SNRdB_vec,ar_data_harq);
hold on;
    semilogy(SNRdB_vec,ar_data_fb);
xlabel('SNR');
ylabel('Avg. rounds');
legend("HARQ","FB");
title_str = sprintf('HARQ vs FB scheme opt. err thr : BLER LDPC %d Rate %.3f max. iter %d max.rounds %d',N,R, max_iter, max_rounds);
title(title_str);
bler_common_str = [res_folder sprintf('/AR_HARQ_vs_FB_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_ada_max_rounds_%d_qm_%d_ma_%d_numF_%d',N,R, max_iter, max_rounds, qam_mod, mod_approx, nFrames)];
filename_BLER_fig = bler_common_str + ".fig";
filename_BLER_png = bler_common_str + ".png";

savefig(filename_BLER_fig);
saveas(f,filename_BLER_png);

% % Terminate parpool
% if (parpool_size ~= 8 && parpool_size ~= 10)
%     delete(gcp('nocreate'));
% end