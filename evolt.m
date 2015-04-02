% echarge  Electron volt (unit of energy)
%
%   eV = evolt
%   [eV,sigma] = evolt
%
%   Returns the electron volt, a conventional unit of energy,
%   in SI units (Joule). sigma is the standard uncertainty (2010 CODATA).

function [eV,sigma] = evolt

% 2010 CODATA value
%Value 	 1.602 176 565 x 10-19 C
% Standard uncertainty 	 0.000 000 035 x 10-19 C
%  Relative standard uncertainty 	  2.2 x 10-8
% Concise form 	 1.602 176 565(35) x 10-19 C
 
eV =    1.602176565e-19;
sigma = 0.000000035e-19;

return
