function [C, V, beta] = compute_BIAWGN_params(snr)

% This script is to compute the channel capacity and dispersion for the
% BI-AWGN channel at the given SNR. The expressions for each quantity can
% be found in research slides dated 01-05-22, "Theory of SDO".

% Input parameters:
%   1) snr: the E_s/(\sigma^2) in dB form.
%
% Output parameters:
%   1) C: the capacity of the BI-AWGN channel at the given SNR (bits)
%   2) V: the dispersion (variance) of the BI-AWGN channel at the given SNR
%   3) beta: the absolute 3rd central moment
%
% Remarks:
%   1) The logarithms are base 2.
%
% Written by Hengjie Yang (hengjie.yang@ucla.edu)   01/07/22
%

P = 10^(snr/10);
zmin = -10;
zmax = 10;

C = integral(@(z) 1/sqrt(2*pi).*exp(-z.^2./2).*(1 - log2(1 + exp(-2*P-2*sqrt(P).*z))), zmin, zmax);
V = integral(@(z) 1/sqrt(2*pi).*exp(-z.^2./2).*(1 - log2(1 + exp(-2*P-2*sqrt(P).*z))).^2, zmin, zmax);
V = V - C^2;
beta = integral(@(z) 1/sqrt(2*pi).*exp(-z.^2./2).*abs((1 - log2(1 + exp(-2*P-2*sqrt(P).*z)) - C).^3), zmin, zmax);




end