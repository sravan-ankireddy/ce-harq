function [ys] = Petrov_CDF(snr, n, gamma)

% The function applies the order-3 Petrov expansion to approximate the CDF
% of the cumulative information density for the BI-AWGN channel. An
% explicit expression of the expansion can be found in the research slides
% dated 01-05-22 in folder "Theory of Sequential Differential
% Optimization".
%
% Input parameters:
%   1) snr: a scalar denoting the BI-AWGN channel SNR in dB.
%   2) n: a scalar (row vector) denoting the number of i.i.d. RVs in the
%       cumulative information density.
%   3) gamma: a row vector (scalar) denoting the \gamma of interest.
%
% Output parameters:
%   1) ys: a row vector denoting Pr(i(X^n; Y^n) <= \gamma).
%
% Remarks:
%   1) This function relies on the following functions:
%       (1) compute_BIAWGN_params.m
%
%   2) The vectorized computation is supported in the following cases:
%       (1) a scalar 'n', and a vector 'gamma';
%       (2) a vector 'n', and a scalar 'gamma';
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/18/22
%


% basic parameters
P = 10^(snr/10); % the BPSK signal power
[C, V] = compute_BIAWGN_params(snr);
zs = (gamma - n.*C)./sqrt(n.*V); % the normalized input value
N = 3; % the book only provided order-3 Petrov expansion
nonctr_moments = zeros(1, N+2); % noncentral moments E[Y^l] up to order N+2
kappas = zeros(1, N+2); % cumulants of (X-C) up to order N+2

% compute noncentral moments E[Y^m]
for m = 2:(N+2)
    nonctr_moments(m) = integral(@(z) 1/sqrt(2*pi).*exp(-z.^2/2).*(1-log2(1+exp(-2*P-2*sqrt(P).*z))-C).^m, -10, 10);
end

% compute cumulants \kappa_m
for m = 1:(N+2)
    [nsol, Sols] = find_solutions(m);
    for ii = 1:nsol
        r = sum(Sols(ii, :));
        temp = 1;
        for ll = 1:m
            temp = temp * 1/factorial(Sols(ii, ll)) * (nonctr_moments(ll)/factorial(ll))^(Sols(ii, ll));
        end
        kappas(m) = kappas(m) + (-1)^(r-1) * factorial(r-1) * temp;
    end
    kappas(m) = factorial(m)*kappas(m);
end


% break zs into non-negative and negative parts
temps = zs.^3./sqrt(n);


a0 = kappas(3) / (6*kappas(2)^(3/2));
a1 = (kappas(4)*kappas(2) - 3*kappas(3)^2) / (24*kappas(2)^3);
a2 = (kappas(5)*kappas(2)^2 - 10*kappas(4)*kappas(3)*kappas(2) + 15*kappas(3)^3) / (120*kappas(2)^(9/2));

% compute ys for zs non-negative
ys_pos = ones(size(zs));
lambdas = a0 + a1.*zs./sqrt(n) + a2.*(zs./sqrt(n)).^2;
mid_seq = temps.*lambdas;
ys_pos = ys_pos - qfunc(zs).*exp(temps.*lambdas);

% a = 0.339;
% b = 5.51;
% mid_seq = -zs.^2/2 + temps.*lambdas;
% ys_pos = ys_pos - vpa(1/sqrt(2*pi)*1./((1-a).*zs + a*sqrt(zs.^2 + b)).*exp(mid_seq));

ys_pos = ys_pos(zs >= 0);

% compute ys for zs negative
ys_neg = zeros(size(zs));
lambdas = a0 + a1.*(zs./sqrt(n)) + a2.*(zs./sqrt(n)).^2;
ys_neg = ys_neg + normcdf(zs).*exp(temps.*lambdas);
ys_neg = ys_neg(zs < 0);

ys = [ys_pos, ys_neg];


end





