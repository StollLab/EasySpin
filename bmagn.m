% bmagn  Bohr magneton 
%
%   mu = bmagn
%   [mu,sigma] = bmagn
%
%   Returns the Bohr magneton in Joule/Tesla.
%   sigma is the standard uncertainty (2006
%   CODATA).

function [mu,sigma] = bmagn

% CODATA 2006 value
%Value 927.400 915 x 10-26 J T-1
%Standard uncertainty 0.000 023 x 10-26 J T-1
%Relative standard uncertainty 2.5 x 10-8
%Concise form 927.400 915(23) x 10-26 J T-1

mu =    9.27400915e-24;
sigma = 0.00000023e-24;
