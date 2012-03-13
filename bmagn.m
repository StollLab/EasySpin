% bmagn  Bohr magneton 
%
%   mu = bmagn
%   [mu,sigma] = bmagn
%
%   Returns the Bohr magneton in Joule/Tesla.
%   sigma is the standard uncertainty (2010
%   CODATA).

function [mu,sigma] = bmagn

% CODATA 2010 value 927.400 968(20) x 10^-26 J T-1

mu =    9.27400968e-24;
sigma = 0.00000020e-24;
