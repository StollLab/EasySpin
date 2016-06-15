% boltzm  Boltzmann constant 
%
%   k = boltzm
%   [k,sigma] = boltzm
%
%   Returns the Boltzmann constant in SI units, joule per kelvin.
%   sigma is the standard uncertainty.

function [k_B,sigma] = boltzm

% 2014 CODATA value
%
% Boltzmann constant
% Concise form 	 1.380 648 52(79) x 10-23 J K-1

k_B =   1.38064852e-23;
sigma = 0.00000079e-23;
