% echarge  Electron volt (unit of energy)
%
%   eV = evolt
%   [eV,sigma] = evolt
%
%   Returns the electron volt, a conventional unit of energy,
%   in SI units (Joule). sigma is the standard uncertainty.

function [eV,sigma] = evolt

% 2014 CODATA value
% Concise form 	 1.602 176 6208(98) x 10-19 C
 
eV =    1.6021766208e-19;
sigma = 0.0000000098e-19;
