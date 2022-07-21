close all; clear all; clc;

global_settings = 1;

unit_test = 0;

% run params
gen_err_thr_harq = 0;
gen_err_thr_fb = 0;
gen_harq_vs_fb = 1;

process_data_harq = 1;
process_data_fb = 1;

% LUT path
err_thr_lut_path = "lut_data/fb_lut_data_LDPC_960_rate_0.904_dec_iter_6_err_thr_0.010_to_0.200_max_rounds_10_qm_0_ma_0_numF_10000.mat";
harq_lut_path = "lut_data/harq_lut_data_LDPC_960_dec_iter_6_rate_0.050_to_0.900_max_rounds_10_qm_0_ma_0_numF_10000.mat";

% generate acomp LUT
acomp_lut_path = "lut_data/acomp_960_ns_100000.mat";
if ~isfile(acomp_lut_path)
   acomp_table = generate_acomp_table;
else
    acomp_table = load(acomp_lut_path);
end

% constant settings
nPRB = 20;
max_iter = 6;
mod_approx = 0;
qam_mod = 0;
max_rounds = 10;
nFrames = 10e3;
min_bler = 1e-4; % URLLC requirement
modulation = 'QPSK';
% HARQ params
targetCodeRate_list = 0.05:0.05:0.9;
actualCodeRate_list = zeros(size(targetCodeRate_list));
SNRdB_low = -20;
SNRdB_high = 10;
SNRdB_step = 0.2;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
num_SNRdB = length(SNRdB_vec);

if (unit_test)
    nFrames = 40e1;
    targetCodeRate_list = 0.05:0.05:0.1;
    actualCodeRate_list = zeros(size(targetCodeRate_list));
end
% Store full error data : 2xtargetCodeRate_listxmax_roundsxnum_SNR
err_data_all_rates = zeros(2,length(targetCodeRate_list),max_rounds,num_SNRdB);

% Store Avg. rounds data : targetCodeRate_listxnum_SNR
ar_data_all_rates = zeros(length(targetCodeRate_list),num_SNRdB);

% Initiate parpool
parpool_size = min(64,feature('numcores'));
if (parpool_size ~= 8 && parpool_size ~= 10)
    parpool(parpool_size);
end

snr_data = SNRdB_vec;

if (gen_err_thr_harq == 1)
    res_folder = sprintf('bler_data/harq/%d',nFrames);

    if ~exist(res_folder,'dir')
        mkdir(res_folder);
    end
    
    res_folder_all = res_folder + "/all";
    if ~exist(res_folder_all,'dir')
        mkdir(res_folder_all);
    end
    err_thr = 0.05;
    err_thr_ada_list = err_thr*ones(size(SNRdB_vec));

    for i_tr = 1:length(targetCodeRate_list)
        
        targetCodeRate = targetCodeRate_list(i_tr);

        % Call the main script
        run_harq;

        % Save the results            
        data_file_name = [res_folder sprintf('/harq_lut_data_LDPC_%d_dec_iter_%d_rate_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d.mat',N, max_iter, targetCodeRate_list(1), targetCodeRate_list(end), max_rounds, qam_mod, mod_approx, nFrames)];

        save(data_file_name,'err_data_all_rates','ar_data_all_rates','snr_data','targetCodeRate_list','actualCodeRate_list');
    end
    %% Plot the cumulative results
    % set the figure properties for BLER plots
    f = figure('Renderer','painters','Position',[1000 400 800 500]);

    leg_str = {};
    for i_tr = 1:length(actualCodeRate_list)
        semilogy(SNRdB_vec,squeeze(err_data_all_rates(2,i_tr,end,:)));
        hold on;
        leg_str{end+1} = sprintf('Rate %.3f',actualCodeRate_list(i_tr));
    end
    xlabel('SNR');
    ylabel('BLER');
    legend(leg_str);
    title_str = sprintf('HARQ scheme : BLER LDPC %d max. rounds %d',N,max_rounds);
    title(title_str);
    bler_common_str = [res_folder sprintf('/HARQ_BLER_LDPC_%d_dec_iter_%d_rate_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d_',N, max_iter, targetCodeRate_list(1), targetCodeRate_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
    filename_BLER_fig = bler_common_str + ".fig";
    filename_BLER_png = bler_common_str + ".png";

    savefig(filename_BLER_fig);
    saveas(f,filename_BLER_png);

    % set the figure properties for AR plots
    f = figure('Renderer','painters','Position',[1000 400 800 500]);

    leg_str = {};
    for i_tr = 1:length(actualCodeRate_list)
        semilogy(SNRdB_vec,ar_data_all_rates(i_tr,:));
        hold on;
        leg_str{end+1} = sprintf('Rate%.3f',actualCodeRate_list(i_tr));
    end
    xlabel('SNR');
    ylabel('Avg. rounds');
    legend(leg_str);
    title_str = sprintf('HARQ scheme : AR LDPC %d max. rounds %d',N,max_rounds);
    title(title_str);
    ar_common_str = [res_folder sprintf('/HARQ_AR_LDPC_%d_dec_iter_%d_rate_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d_',N, max_iter, targetCodeRate_list(1), targetCodeRate_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
    filename_AR_fig = ar_common_str + ".fig";
    filename_AR_png = ar_common_str + ".png";

    savefig(filename_AR_fig);
    saveas(f,filename_AR_png);
end

%% Feedback LUT
if (gen_err_thr_fb == 1)

    % FB params
    targetCodeRate = 0.8;
    res_folder = sprintf('bler_data/fb/%d',nFrames);
    err_thr_list = 0.01:0.01:0.2;
    SNRdB_low = -6;
    SNRdB_high =-1;
    SNRdB_step = 0.2;
    SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
    num_SNRdB = length(SNRdB_vec);
    snr_data = SNRdB_vec;

    if (unit_test)
        err_thr_list = 0.01:0.01:0.02;
    end
            
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
    title_str = sprintf('FB scheme : BLER LDPC %d Rate %.3f max.rounds %d',N,R, max_rounds);
    title(title_str);
    bler_common_str = [res_folder sprintf('/FB_BLER_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d_',N,R, max_iter, err_thr_list(1), err_thr_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
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
    title_str = sprintf('FB scheme : AR LDPC %d Rate %.3f max. rounds %d',N,R, max_rounds);
    title(title_str);
    ar_common_str = [res_folder sprintf('/FB_AR_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d_',N,R, max_iter, err_thr_list(1), err_thr_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
    filename_AR_fig = ar_common_str + ".fig";
    filename_AR_png = ar_common_str + ".png";

    savefig(filename_AR_fig);
    saveas(f,filename_AR_png);
end

%% HARQ vs FB with opt. err_thr
if (gen_harq_vs_fb == 1)
    
    % HARQ vd FB params
    targetCodeRate = 0.9;
    SNRdB_low = -5;
    SNRdB_high = 1;
    SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
    num_SNRdB = length(SNRdB_vec);
    % use either estimate or opt values from LUT
    err_thr_ada_scheme = "est";

    % generate estimates for err_thr
    if (err_thr_ada_scheme == "est")
        data = load(harq_lut_path);
        err_thr_ada_list_est = zeros(max_rounds,length(SNRdB_vec));
        for i_rr = 1:size(err_thr_ada_list_est,1)
            remRounds = i_rr;
            for i_ada = 1:size(SNRdB_vec,2)
                err_thr_ada_list_est(i_rr,i_ada) =  err_thr_select(data,acomp_table,targetCodeRate,SNRdB_vec(i_ada),remRounds,min_bler);
            end
        end
    end

    if (err_thr_ada_scheme == "est")
        res_folder = sprintf('bler_data/harq_vs_fb_est/%d',nFrames);

        if ~exist(res_folder,'dir')
            mkdir(res_folder);
        end
    else
        res_folder = sprintf('bler_data/harq_vs_fb/%d',nFrames);

        if ~exist(res_folder,'dir')
            mkdir(res_folder);
        end
    end

    % look up opt err_thr per SNR
    data = load(err_thr_lut_path);
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

    % HARQ
    % Store full error data : 2xmax_roundsxnum_SNR
    err_data_harq = zeros(2,max_rounds,num_SNRdB);

    % call the main script
    process_data_harq = 0;
    run_harq;   
    err_data_harq(1,:,:) = BER_vec_pr_HARQ;
    err_data_harq(2,:,:) = BLER_vec_pr_HARQ;
    ar_data_harq = Avg_rounds_HARQ;
    
    % FB
    % Store full error data : 2xmax_roundsxnum_SNR
    err_data_fb = zeros(2,max_rounds,num_SNRdB);

    % call the main script
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
    title_str = sprintf('HARQ vs FB scheme opt. err thr : BLER LDPC %d Rate %.3f max.rounds %d',N,R, max_rounds);
    title(title_str);
    bler_common_str = [res_folder sprintf('/BLER_HARQ_vs_FB_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_ada_max_rounds_%d_qm_%d_ma_%d_numF_%d_',N,R, max_iter, max_rounds, qam_mod, mod_approx, nFrames)];
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
    title_str = sprintf('HARQ vs FB scheme opt. err thr : BLER LDPC %d Rate %.3f max.rounds %d',N,R, max_rounds);
    title(title_str);
    bler_common_str = [res_folder sprintf('/AR_HARQ_vs_FB_LDPC_%d_rate_%.3f_dec_iter_%d_err_thr_ada_max_rounds_%d_qm_%d_ma_%d_numF_%d_',N,R, max_iter, max_rounds, qam_mod, mod_approx, nFrames)];
    filename_BLER_fig = bler_common_str + ".fig";
    filename_BLER_png = bler_common_str + ".png";

    savefig(filename_BLER_fig);
    saveas(f,filename_BLER_png);
end


% Terminate parpool
if (parpool_size ~= 8 && parpool_size ~= 10)
    delete(gcp('nocreate'));
end


