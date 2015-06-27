% faraday  Faraday constant 
%
%   F = faraday
%   [F,sigma] = faraday
%
%   Returns the Faraday constant in SI units, coulomb per mole.
%   sigma is the standard uncertainty.

function [F,sigma] = faraday

% 2014 CODATA value
% Concise form 	 96 485.332 89(59) C mol-1 
 
F = 96485.33289;
sigma = 0.00059;
