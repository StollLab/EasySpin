% echarge  Elementary charge 
%
%   e = echarge
%   [e,sigma] = echarge
%
%   Returns the elementary charge in SI units,
%   coulomb. sigma is the standard uncertainty.

function [e,sigma] = echarge

% 2014 CODATA value
% Concise form 	 1.602 176 6208(98) x 10-19 C
 
e =     1.6021766208e-19;
sigma = 0.0000000098e-19;
