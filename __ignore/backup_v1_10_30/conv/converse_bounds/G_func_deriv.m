function [ys] = G_func_deriv(snr, ns, gamma)

% This function computes the first derivative of the order-3 Petrov expansion.
% Since Petrov expansion has two different expressions depending on the
% sign of the input argument, the first derivative also has two different
% expressions depending on the sign of the input argument. Explicit
% derivative expressions can be found in the research slides dated
% 01-05-22 in the folder "Theory of Sequential Differential Optimization".
%
% Input parameters:
%   1) snr: a scalar denoting the BI-AWGN channel SNR in dB.
%   2) ns: a row vector denoting the number of i.i.d. RVs in the
%       cumulative information density. Must be a monotonically increasing
%       sequence.
%   3) gamma: a scalar denoting the \gamma of interest.
%
% Output parameters:
%   1) ys: a row vector denoting the first derivative w.r.t. n.
%
% Remarks:
%   1) This function relies on the following functions:
%       (1) compute_BIAWGN_params.m
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/20/22
%


% basic parameters
P = 10^(snr/10); % the BPSK signal power
[C, V] = compute_BIAWGN_params(snr);
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


seq = (gamma - ns.*C)./sqrt(ns.*V); % the normalized sequence
seq_pos = seq(seq >= 0);
seq_neg = seq(seq < 0);
ns_large = ns(length(seq_pos)+1:end);
ns_small = ns(1:length(seq_pos));

a0 = kappas(3) / (6*kappas(2)^(3/2));
a1 = (kappas(4)*kappas(2) - 3*kappas(3)^2) / (24*kappas(2)^3);
a2 = (kappas(5)*kappas(2)^2 - 10*kappas(4)*kappas(3)*kappas(2) + 15*kappas(3)^3) / (120*kappas(2)^(9/2));


% compute ys associated with n small
ys_small = zeros(size(ns_small));
lambda_small = a0 + a1.*(seq_pos./sqrt(ns_small)) + a2.*(seq_pos./sqrt(ns_small)).^2;
lambda_deriv_small = -(gamma./(ns_small.^2*sqrt(V))).*(a1 + 2*a2.*(gamma - ns_small.*C)./(ns_small.*sqrt(V)));
ys_small = ys_small + (gamma + ns_small.*C)./(2*ns_small.^(3/2).*sqrt(V)).*normpdf(seq_pos).*...
    exp((gamma - ns_small.*C).^3./(ns_small.^2.*V^(3/2)).*lambda_small) + ...
    qfunc(seq_pos).*exp((gamma - ns_small.*C).^3./(ns_small.^2.*V^(3/2)).*lambda_small).*...
    (-(gamma - ns_small.*C).^2.*(ns_small*C+2*gamma)./(ns_small.^3.*V^(3/2)).*lambda_small +...
    (gamma - ns_small.*C).^3./(ns_small.^2.*V^(3/2)).*lambda_deriv_small );


% compute ys associated with n large
ys_large = zeros(size(ns_large));
lambda_large = a0 + a1.*(seq_neg./sqrt(ns_large)) + a2.*(seq_neg./sqrt(ns_large)).^2;
lambda_deriv_large = -(gamma./(ns_large.^2*sqrt(V))).*(a1 + 2*a2.*(gamma - ns_large.*C)./(ns_large.*sqrt(V)));
ys_large = ys_large + (gamma + ns_large.*C)./(2*ns_large.^(3/2).*sqrt(V)).*normpdf(seq_neg).*...
    exp((gamma - ns_large.*C).^3./(ns_large.^2.*V^(3/2)).*lambda_large) - ...
    normcdf(seq_neg).*exp((gamma - ns_large.*C).^3./(ns_large.^2.*V^(3/2)).*lambda_large).*...
    (-(gamma - ns_large.*C).^2.*(ns_large*C+2*gamma)./(ns_large.^3.*V^(3/2)).*lambda_large +...
    (gamma - ns_large.*C).^3./(ns_large.^2.*V^(3/2)).*lambda_deriv_large );


% merge results together
ys = [ys_small, ys_large];



end





