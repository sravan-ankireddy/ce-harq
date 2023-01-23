function [ys] = F_func_deriv_def(snr, ns, s, gamma)

% This script is to use the definition of derivative to compute f function.
% The purpose of this script is to cross-verify the original first
% derivative script. The idea is to use F(x+\delta) - F(x)/\delta to
% approximate the derivative at x.
%
% Input parameters:
%   1) snr: a scalar denoting the BI-AWGN channel SNR in dB;
%   2) ns: a row vector indicating the n's of interest;
%   3) s: a scalar denoting the order of Edgeworth expansion
%   4) gamma: a scalar denoting the cumu. info. density threshold
%
% Output parameters:
%   1) ys: a row vector denoting dF/dn, where F = 1 - Edgeworth_CDF.
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/21/22
%


delta = 1e-6;
F_vals = 1 - Edgeworth_CDF(snr, ns, s, gamma);
ns_shifted = ns + delta;
F_val_shifted = 1 - Edgeworth_CDF(snr, ns_shifted, s, gamma);
ys = (F_val_shifted - F_vals)./(ns_shifted - ns);



end