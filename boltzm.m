% boltzm  Boltzmann constant 
%
%   k = boltzm
%   [k,sigma] = boltzm
%
%   Returns the Boltzmann constant in Joule/Kelvin.
%   sigma is the standard uncertainty (2010 CODATA).

function [k_b,sigma] = boltzm

% 2010 CODATA value
%
% Boltzmann constant
%Value 	 1.380 6488 x 10-23 J K-1
%Standard uncertainty 	 0.000 0013 x 10-23 J K-1
%  Relative standard uncertainty 	  9.1 x 10-7
% Concise form 	 1.380 6488(13) x 10-23 J K-1

k_b =   1.3806488e-23;
sigma = 0.0000013e-23;

return
