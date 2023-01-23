function [ys] = G_over_g_ratio(snr, ms, ns, gamma)

% This function is to evaluate G(m)/g(n), where G is the tail probability
% approximated by Petrov expansion and g = dG/dn.
%
% Input parameters:
%   1) snr: a scalar denoting the BI-AWGN channel SNR in dB.
%   2) ms: a length-N row vector denoting the input argument of G(m); m < \gamma/C
%   3) ns: a length-N row vector denoting the input argument of g(n); n < \gamma/C
%   4) gamma: a scalar denoting the \gamma of interest.
%
% Output parameters:
%   1) ys: a row vector denoting G(m)/g(n).
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/21/22
%


% basic parameters
P = 10^(snr/10);
[C, V] = compute_BIAWGN_params(snr);
N = 3; % the book only provided order-3 Petrov expansion
nonctr_moments = zeros(1, N+2); % noncentral moments E[Y^l] up to order N+2
kappas = zeros(1, N+2); % cumulants of (X-C) up to order N+2

% check whether the input meets the requirement
if max(ms) > gamma/C
    error('The input m exceeds \gamma/C!');
elseif max(ns) > gamma/C
    error('The input n exceeds \gamma/C!');
end

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

a0 = kappas(3) / (6*kappas(2)^(3/2));
a1 = (kappas(4)*kappas(2) - 3*kappas(3)^2) / (24*kappas(2)^3);
a2 = (kappas(5)*kappas(2)^2 - 10*kappas(4)*kappas(3)*kappas(2) + 15*kappas(3)^3) / (120*kappas(2)^(9/2));


x1s = (gamma - ms.*C)./(ms.*sqrt(V));
lambda1s = a0 + a1.*x1s + a2.*x1s.^2;
z1s = (gamma - ms.*C).^3./(ms.^2.*V.^(3/2)).*lambda1s;



xs = (gamma  - ns.*C)./(ns.*sqrt(V));
lambdas = a0 + a1.*xs + a2.*xs.^2;
lambda_derivs = -gamma./(ns.^2.*sqrt(V)).*(a1 + 2*a2.*xs);
zs = (gamma - ns.*C).^3./(ns.^2.*V.^(3/2)).*lambdas;

diff_zs = z1s - zs;
ys = qfunc(x1s.*sqrt(ms)).*exp(diff_zs)./((gamma+ns.*C)./(2*ns.^(3/2).*sqrt(V)).*normpdf(xs.*sqrt(ns)) + ...
    qfunc(xs.*sqrt(ns)).*(-(gamma - ns.*C).^2.*(ns.*C+2*gamma)./(ns.^3.*V^(3/2)).*lambdas + ...
    (gamma - ns.*C).^3./(ns.^2.*V^(3/2)).*lambda_derivs ) );




end




