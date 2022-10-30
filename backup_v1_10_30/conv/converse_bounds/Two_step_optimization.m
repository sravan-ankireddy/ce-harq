% This script is to implement the two-step optimization over (\gamma,
% n_1^m). First, we search over \gamma. Next, for a given \gamma, we use
% SDO to find conditionally optimal decoding times n_1^m and compute its
% average blocklength. The \gamma that yields the minimum average
% blocklength is the optimal \gamma. We would like to compare this
% achivability bound with Polyanskiy's optimized VLSF bound:
%   E[\tau] <= (\log(M-1) - \log\epsilon + 1) / C.
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/30/22
%


clear;
clc;


tic
% basic parameters
snr = 0.2; % E_s/(sigma^2) in dB
ks = [1:20]; % information lengths
ms = [1:8]; % number of decoding opportunities
s = 5; % order of Edgeworth expansion
epsilon = 1e-3; % target error probability


%% Plot the optimal average blockelength vs. delta and see if it is convex

% deltas = (1:99)/100;
% vals = zeros(size(deltas));
% parfor ii = 1:length(deltas)
%     disp(['Current progress: ii = ', num2str(ii)]);
%     vals(ii) = find_ave_bl(snr, epsilon, 20, s, 3, deltas(ii));
% end
% plot(deltas, vals, 'r-');
% grid on
% xlabel('$\delta$','Interpreter','latex');
% ylabel('Average blocklength');



%%
Optimal_deltas = zeros(length(ms), length(ks));
Optimal_dec_times = cell(length(ms), length(ks));
Optimal_ave_blocklengths = zeros(length(ms), length(ks));



for iter = 1:length(ms)
    m = ms(iter);
%     m = 3;
    parfor ii = 1:length(ks)
        k = ks(ii);
        disp(['Current progress: iter = ', num2str(iter), ' out of ',...
            num2str(length(ms)),', ii = ', num2str(ii), ' out of ', num2str(length(ks))]);
        
        tolerance = 1e-3;
        upper_delta = 1 - tolerance;
        lower_delta = tolerance;

        opt_delta = fminbnd(@(z) find_ave_bl(snr, epsilon, k, s, m, z), tolerance, 1-tolerance);
        
%         while abs(upper_delta - lower_delta) > tolerance
%             mid_delta = (lower_delta + upper_delta) / 2;
%             Ave_bl_upper = find_ave_bl(snr, epsilon, k, s, m, upper_delta);
%             Ave_bl_lower = find_ave_bl(snr, epsilon, k, s, m, lower_delta);
%             if Ave_bl_upper <= Ave_bl_lower
%                 lower_delta = mid_delta;
%             else
%                 upper_delta = mid_delta;
%             end
%         end
%         opt_delta = upper_delta;
        [opt_bl, opt_Ns] = find_ave_bl(snr, epsilon, k, s, m, opt_delta);
        Optimal_deltas(iter, ii) = opt_delta;
        Optimal_dec_times{iter, ii} = opt_Ns;
        Optimal_ave_blocklengths(iter, ii) = opt_bl;
    end
end



%% save dataset


path = './sim_data/';
timestamp = datestr(now, 'mmddyy_HHMMSS_');
file_name = [path, timestamp,'dataset_two_step_optimization_snr_',num2str(snr),'_s_',num2str(s),'_real_k_100','.mat'];
% save(file_name, 'snr', 'epsilon', 'ks', 'ms', 's', 'Optimal_ave_blocklengths', 'Optimal_dec_times', 'Optimal_deltas');






toc


function [ave_bl, Ns] = find_ave_bl(snr, epsilon, k, s, m, delta)


C = compute_BIAWGN_params(snr);
gamma = log2((2^k-1)/(delta*epsilon));
n_m = fzero(@(z) Edgeworth_CDF(snr, z, s, gamma) - epsilon*(1 - delta), max(300, k/C));
[Ns, Ps, ~] = bisec_search_SDO(n_m, snr, s, gamma, m);

ave_bl = Ns(1);
for kk = 1:(m-1)
    ave_bl = ave_bl + (Ns(kk+1) - Ns(kk))*Ps(kk);
end

end












