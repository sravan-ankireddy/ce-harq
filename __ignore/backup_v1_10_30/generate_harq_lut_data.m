close all; clear all; clc;

% Figure settings
startup;

global_settings = 1;

unit_test = 0;

% run params
generate_harq_lut = 1;


% constant settings
% FIX ME
ncb = 1;
Nref = 25344;
max_iter = 6;
mod_approx = 0;
qam_mod = 0;
max_rounds = 10;
nFrames = 10e3;
min_bler = 1e-4; % URLLC requirement 1e-5
modulation = 'QPSK';

% use nrTBS to get K,R
nPRB = 20;
nlayers = 1;
NREPerPRB = 12*4; % For URLLC, 2-7 is the typical choice

% CodeLen based on PRB settings
N = nPRB*NREPerPRB;

% HARQ params
targetCodeRate_list = 0.05:0.05:0.9;
actualCodeRate_list = zeros(size(targetCodeRate_list));
SNRdB_low = -20;
SNRdB_high = 0;
SNRdB_step = 0.2;
SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
num_SNRdB = length(SNRdB_vec);

if (unit_test)
    SNRdB_low = -20;
    SNRdB_high = 10;
    SNRdB_step = 0.2;
    SNRdB_vec = SNRdB_low:SNRdB_step:SNRdB_high;
    num_SNRdB = length(SNRdB_vec);
    nFrames = 10e2;
    targetCodeRate_list = 0.9;%0.05:0.05:0.1;
    actualCodeRate_list = zeros(size(targetCodeRate_list));
end

% Choose the combining scheme 
combining_scheme = "CC";
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

res_folder_harq = [res_folder_prefix sprintf('/%d/harq/%d',N,nFrames)];

% Store full error data : 2xtargetCodeRate_listxmax_roundsxnum_SNR
err_data_all_rates = zeros(2,length(targetCodeRate_list),max_rounds,num_SNRdB);

% Store Avg. rounds data : targetCodeRate_listxnum_SNR
ar_data_all_rates = zeros(length(targetCodeRate_list),num_SNRdB);

% Initiate parpool
parpool_size = min(64,feature('numcores'));
pp = gcp('nocreate'); % If no pool, do not create new one.
if (parpool_size ~= 8 && parpool_size ~= 10 && isempty(pp))
    parpool(parpool_size);
end

snr_data = SNRdB_vec;

if (generate_harq_lut == 1)
    res_folder = res_folder_harq;

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
        process_data_harq = 1;
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
    title_str = sprintf('HARQ %s scheme : BLER LDPC %d max. iter %d max. rounds %d', combining_scheme, N,max_iter,max_rounds);
    title(title_str);
    bler_common_str = [res_folder sprintf('/HARQ_BLER_LDPC_%d_dec_iter_%d_rate_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d',N, max_iter, targetCodeRate_list(1), targetCodeRate_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
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
    title_str = sprintf('HARQ %s scheme : AR LDPC %d max. iter %d max. rounds %d',combining_scheme,N,max_iter,max_rounds);
    title(title_str);
    ar_common_str = [res_folder sprintf('/HARQ_AR_LDPC_%d_dec_iter_%d_rate_%.3f_to_%.3f_max_rounds_%d_qm_%d_ma_%d_numF_%d',N, max_iter, targetCodeRate_list(1), targetCodeRate_list(end), max_rounds, qam_mod, mod_approx, nFrames)];
    filename_AR_fig = ar_common_str + ".fig";
    filename_AR_png = ar_common_str + ".png";

    savefig(filename_AR_fig);
    saveas(f,filename_AR_png);
end

% % Terminate parpool
% if (parpool_size ~= 8 && parpool_size ~= 10)
%     delete(gcp('nocreate'));
% end


