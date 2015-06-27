% bmagn  Bohr magneton 
%
%   mu = bmagn
%   [mu,sigma] = bmagn
%
%   Returns the Bohr magneton in joule per tesla.
%   sigma is the standard uncertainty.

function [mu,sigma] = bmagn

% CODATA 2014 value 927.400 9994(57) x 10^-26 J T-1

mu =    9.274009994e-24;
sigma = 0.000000057e-24;
