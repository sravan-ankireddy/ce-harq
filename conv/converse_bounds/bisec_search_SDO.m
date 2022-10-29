function [Ns, Ps, gap] = bisec_search_SDO(n_m, snr, s, gamma, m)

% The function applies bisection search to identify n_1^*, which
% subsequently yields n_2^*, ..., n_m^* with new SDO procedure. 
%
% Input parameters:
%   1) n_m: a scalar denoting real-valued n_m
%   2) snr: a scalar denoting BI-AWGN channel SNR in dB
%   3) s: a scalar denoting the order of Edgeworth expansion
%   4) gamma: a scalar denoting the cumu. info. density threshold
%   5) m: a scalar denoting the number of decoding times
%
% Output parameters
%   1) the length of Ns determines the number m' of distinct integers after
%       rounding.
%   2) Ps correspond to the integer solution.
%
% 01-22-22 update: This function is a copy from
% "plot_decoding_times_evolution".
%

F_func = @(n) 1 - Edgeworth_CDF(snr, n, s, gamma);
f_func = @(n) F_func_deriv_def(snr, n, s, gamma); % dF(n)/dn

G_func = @(n) 1 - Petrov_CDF(snr, n, gamma);
g_func = @(n) G_func_deriv(snr, n, gamma); % dG(n)/dn

Ns = zeros(1,m); % a row vector consisting of optimal real-valued n_1, n_2, ..., n_m



% identify x_0
x_mid = fzero(@(z) F_func(z) - 1/3, n_m); % identify F(x_mid) = 1/3.
x_0 = fzero(@(z) F_func(z) - G_func(z), x_mid); % use x_mid to identify F(x^*) = G(x^*).


lower = 0.5; % this guarantees the rounding is at least 1.
upper = ceil(n_m) - m + 0.5;
tolerance = 1e-5;

if lower > upper
    error('m is too large that makes the last zero-crossing value fail!');
end

while abs(upper - lower) > tolerance
    mid = (lower + upper) / 2;
    a = 0; % the (k-1)th decoding time
    b = mid; % the kth decoding time
    c = -1; % the (k+1)th decoding time
    prev_lambda_r = 0; % \lambda_{k-1}
    lambda_r = -1; % \lambda_k
    for kk = 1:(m-1)
        if kk == 1
            if b >= x_0 && ~isnan(b)
                c = b + max(1, F_func(b)/f_func(b));
                lambda_r = max(1 - F_func(b)/f_func(b), 0);
            elseif ~isnan(b)
                c = b + max(1, G_over_g_ratio(snr, b, b, gamma));
                lambda_r = max(1 - G_over_g_ratio(snr, b, b, gamma), 0);
            else
                c = b; % special processing for NaN
                lambda_r = -1;
            end
        else % a possible issue is to decide whether a < x_0 and b < x_0
            if b >= x_0 && ~isnan(b)
%                 c = b + max(1, (F_func(b) - F_func(a)) / f_func(b) - prev_lambda_r*f_func(a)/f_func(b) );
%                 lambda_r = max(prev_lambda_r*f_func(a)/f_func(b) + 1 - F_func(b)/f_func(b) + F_func(a)/f_func(b), 0 );
                if prev_lambda_r + f_func(b)/f_func(a) - F_func(b)/f_func(a) + F_func(a)/f_func(a) < 0
                    c = b + (F_func(b) - F_func(a)) / f_func(b) - prev_lambda_r*f_func(a)/f_func(b);
                    lambda_r = 0;
                else
                    c = b + 1;
                    lambda_r = prev_lambda_r*f_func(a)/f_func(b) + 1 - F_func(b)/f_func(b) + F_func(a)/f_func(b);
                end
%             elseif b >= x_0 && ~isnan(b) && a < x_0
%                 if prev_lambda_r + f_func(b)/g_func(a) - F_func(b)/g_func(a) + F_func(a)/g_func(a) < 0
%                     c = b + (F_func(b) - G_func(a)) / f_func(b) - prev_lambda_r*g_func(a)/f_func(b);
%                     lambda_r = 0;
%                 else
%                     c = b + 1;
%                     lambda_r = prev_lambda_r*g_func(a)/f_func(b) + 1 - F_func(b)/f_func(b) + G_func(a)/f_func(b);
%                 end
            elseif ~isnan(b) % b < x_0
%                 c = b + max(1, G_over_g_ratio(snr, b, b, gamma) - G_over_g_ratio(snr, a, b, gamma) - prev_lambda_r * g_ratio(snr, a, b, gamma) );
%                 lambda_r = max(prev_lambda_r * g_ratio(snr, a, b, gamma) + 1 - G_over_g_ratio(snr, b, b, gamma) + G_over_g_ratio(snr, a, b, gamma), 0);
                if prev_lambda_r + g_ratio(snr, b, a, gamma) - G_over_g_ratio(snr, b, a, gamma) + G_over_g_ratio(snr, a, a, gamma) < 0
%                   c = b + (G_func(b) - G_func(a) - prev_lambda) / g_func(b);
                    c = b + G_over_g_ratio(snr, b, b, gamma) - G_over_g_ratio(snr, a, b, gamma) - prev_lambda_r * g_ratio(snr, a, b, gamma);
                    lambda_r = 0;
                else
                    c = b + 1;
                    lambda_r = prev_lambda_r * g_ratio(snr, a, b, gamma) + 1 - G_over_g_ratio(snr, b, b, gamma) + G_over_g_ratio(snr, a, b, gamma);
                end
            else % b == NaN
                c = b; % special processing for NaN
                lambda_r = -1;
            end
        end
        a = b;
        b = c;
        c = -1;
        prev_lambda_r = lambda_r;
        lambda_r = -1;
    end
    if b < n_m
        lower = mid;
    else
        upper = mid;
    end
end

% disp(['Bisection search completed!', ' m = ',num2str(m)]);
Ns(1) = upper;
a = 0; % the (k-1)th 
b = Ns(1); % the kth
c = -1; % the (k+1)th
prev_lambda_r = 0; % \lambda_{k-1}
lambda_r = -1; % \lambda_k
for kk = 1:(m-1)
    if kk == 1
        if b >= x_0 && ~isnan(b)
            c = b + max(1, F_func(b)/f_func(b));
            lambda_r = max(1 - F_func(b)/f_func(b), 0);
        elseif ~isnan(b)
            c = b + max(1, G_over_g_ratio(snr, b, b, gamma));
            lambda_r = max(1 - G_over_g_ratio(snr, b, b, gamma), 0);
        else
            c = b; % special processing for NaN
            lambda_r = -1;
        end
    else
        if b >= x_0 && ~isnan(b)
%             c = b + max(1, (F_func(b) - F_func(a)) / f_func(b) -  prev_lambda_r*f_func(a)/f_func(b) );
%             lambda_r = max(prev_lambda_r*f_func(a)/f_func(b) + 1 - F_func(b)/f_func(b) + F_func(a)/f_func(b), 0 );
            if prev_lambda_r + f_func(b)/f_func(a) - F_func(b)/f_func(a) + F_func(a)/f_func(a) < 0
                c = b + (F_func(b) - F_func(a)) / f_func(b) -  prev_lambda_r*f_func(a)/f_func(b);
                lambda_r = 0;
            else
                c = b + 1;
                lambda_r = prev_lambda_r*f_func(a)/f_func(b) + 1 - F_func(b)/f_func(b) + F_func(a)/f_func(b);
            end
%         elseif b >= x_0 && ~isnan(b) && a < x_0
%             if prev_lambda_r + f_func(b)/g_func(a) - F_func(b)/g_func(a) + F_func(a)/g_func(a) < 0
%                 c = b + (F_func(b) - G_func(a)) / f_func(b) - prev_lambda_r*g_func(a)/f_func(b);
%                 lambda_r = 0;
%             else
%                 c = b + 1;
%                 lambda_r = prev_lambda_r*g_func(a)/f_func(b) + 1 - F_func(b)/f_func(b) + G_func(a)/f_func(b);
%             end
        elseif ~isnan(b)
%             c = b + max(1, G_over_g_ratio(snr, b, b, gamma) - G_over_g_ratio(snr, a, b, gamma) - prev_lambda_r * g_ratio(snr, a, b, gamma) );
%             lambda_r = max(prev_lambda_r * g_ratio(snr, a, b, gamma) + 1 - G_over_g_ratio(snr, b, b, gamma) + G_over_g_ratio(snr, a, b, gamma), 0);
            if prev_lambda_r + g_ratio(snr, b, a, gamma) - G_over_g_ratio(snr, b, a, gamma) + G_over_g_ratio(snr, a, a, gamma) < 0
%                     c = b + (G_func(b) - G_func(a) - prev_lambda) / g_func(b);
                c = b + G_over_g_ratio(snr, b, b, gamma) - G_over_g_ratio(snr, a, b, gamma) - prev_lambda_r * g_ratio(snr, a, b, gamma);
                lambda_r = 0;
            else
                c = b + 1;
                lambda_r = prev_lambda_r * g_ratio(snr, a, b, gamma) + 1 - G_over_g_ratio(snr, b, b, gamma) + G_over_g_ratio(snr, a, b, gamma);
            end
        else
            c = b; % special processing for NaN
            lambda_r = -1;
        end
    end
    Ns(kk+1) = c;
    a = b;
    b = c;
    c = -1;
    prev_lambda_r = lambda_r;
    lambda_r = -1;
end

% compute Ps
for ii = 1:m
    if Ns(ii) < x_0
        Ps(ii) = 1 - G_func(Ns(ii));
    else
        Ps(ii) = 1 - F_func(Ns(ii));
    end
end

Ns_shifted = Ns(2:end);
gap = min(Ns_shifted - Ns(1:end-1)); % gap < 1 should never happen now.


end
