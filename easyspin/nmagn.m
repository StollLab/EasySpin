% nmagn  Nuclear magneton 
%
%   mu = nmagn
%   [mu,sigma] = nmagn
%
%   Returns the nuclear magneton in SI units, joule per tesla.
%   sigma is the standard uncertainty.

function [mu_n,sigma] = nmagn

% 2014 CODATA recommended values
% Concise form 	 5.050 783 699(31) x 10-27 J T-1

mu_n =  5.050783699e-27;
sigma = 0.000000031e-27;
