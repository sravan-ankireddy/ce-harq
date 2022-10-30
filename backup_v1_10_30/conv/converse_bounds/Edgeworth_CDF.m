function [ys] = Edgeworth_CDF(snr, n, N, gamma)

% The function applies an order-N Edgeworth expansion to approximate the
% CDF of the cumulative information density for the BI-AWGN channel.
% An explicit expression of the expansion can be found in research slides
% dated 01-05-22 in folder "Theory of Sequential Differential
% Optimization".
%
% Input parameters:
%   1) snr: a scalar denoting the BI-AWGN channel SNR in dB.
%   2) n: a scalar (row vector) denoting the number of i.i.d. RVs in the cumulative
%       info. density.
%   3) N: a scalar denoting the order of Edgeworth expansion.
%   4) gamma: a row vector (scalar) denoting the \gamma of interest.
%
% Output parameters:
%   4) ys: a row vector denoting Pr(i(X^n; Y^n) <= \gamma).
%
% Remarks:
%   1) This function relies on the following functions:
%       (1) compute_BIAWGN_params.m
%       (2) find_solutions.m
%   2) The vectorized computation is supported in the following cases:
%       (1) a scalar 'n', and a vector 'gamma';
%       (2) a vector 'n', and a scalar 'gamma';
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/12/22
%



% basic parameters
P = 10^(snr/10); % the BPSK signal power
[C, V] = compute_BIAWGN_params(snr);
zs = (gamma - n.*C)./sqrt(n.*V); % the normalized input value
nonctr_moments = zeros(1, N+2); % noncentral moments E[Y^l] up to order N+2
kappas = zeros(1, N+2); % cumulants of (X-C)/\sqrt(V) up to order N+2

ys = normcdf(zs);

% compute noncentral moments E[Y^m]
for m = 2:(N+2)
    nonctr_moments(m) = integral(@(z) 1/sqrt(2*pi).*exp(-z.^2/2).*((1-log2(1+exp(-2*P-2*sqrt(P).*z))-C)./sqrt(V)).^m, -10, 10);
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


% compute each p_j polynomial

He_inner_func = @(xs, k, j) factorial(j)*(-1)^k.*xs.^(j-2*k)./(factorial(k)*factorial(j-2*k)*2^k);

for jj = 1:N
    [nsol, Sols] = find_solutions(jj);
    p_j_seq = zeros(size(zs));
    for ii = 1:nsol
        r = sum(Sols(ii, :));
        temp = 1;
        for m = 1:jj
            temp = temp * 1/factorial(Sols(ii, m)) * (kappas(m+2)/factorial(m+2))^(Sols(ii, m));
        end

        % initialize the computation of Hermite function of order jj+2*r-1
        He_seq = zeros(size(zs));
        for kk = 0:floor((jj+2*r-1)/2)
            He_seq = He_seq + He_inner_func(zs, kk, jj+2*r-1); % pile up row vectors
        end

        % update the p_j function
        p_j_seq = p_j_seq - He_seq * temp;
    end

    % finally, update ys sequence
    ys = ys + n.^(-jj/2).*normpdf(zs).*p_j_seq;
end

























end