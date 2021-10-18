% nmagn  Nuclear magneton 
%
%   mu = nmagn
%   [mu,sigma] = nmagn
%
%   Returns the nuclear magneton in SI units, joule per tesla.
%   sigma is the standard uncertainty.

function [mu_n,sigma] = nmagn

% 2018 CODATA recommended values

mu_n =  5.0507837461e-27;
sigma = 0.0000000015e-27;
