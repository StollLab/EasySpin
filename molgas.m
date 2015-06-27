% molgas  molar gas constant 
%
%   R = molgas
%   [R,sigma] = molgas
%
%   Returns the molar gas constant in SI units, joule per mole per kelvin.
%   sigma is the standard uncertainty.

function [R,sigma] = molgas

% 2014 CODATA value
% Concise form 	 8.314 4598(48) J mol-1 K-1 

R =     8.3144598;
sigma = 0.0000048;
