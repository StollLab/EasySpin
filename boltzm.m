% boltzm  Boltzmann constant 
%
%   k = boltzm
%   [k,sigma] = boltzm
%
%   Returns the Boltzmann constant in Joule/Kelvin.
%   sigma is the standard uncertainty (2006 CODATA).

function [k_b,sigma] = boltzm

% 2006 CODATA value
%
% Boltzmann constant
% Value 1.380 6504 x 10-23 J K-1
% Standard uncertainty 0.000 0024 x 10-23 J K-1
% Relative standard uncertainty 1.7 x 10-6
% Concise form 1.380 6504(24) x 10-23 J K-1

k_b =   1.3806504e-23;
sigma = 0.0000024e-23;

return
