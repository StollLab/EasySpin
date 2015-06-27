% planck  Planck constant 
%
%   h = planck
%   [h,sigma] = planck
%
%   Returns the Planck constant in SI units, joule times second.
%   sigma is the standard uncertainty.

function [h,sigma] = planck

% 2014 CODATA recommended value
% Concise form 	 6.626 070 040(81) x 10-34 J s

h =     6.626070040e-34;
sigma = 0.000000081e-34;
