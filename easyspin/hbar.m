% hbar    Planck constant over 2 pi
%
%   value = hbar
%   [value,sigma] = hbar
%
%   Returns the Planck constant over 2 pi, h/(2*pi), in joule times second.
%   sigma is the standard uncertainty.

function [value,sigma] = hbar

% 2014 CODATA recommended value
% Concise form 	 1.054 571 800(13) x 10-34 J s
 
value = 1.054571800e-34;
sigma = 0.000000013e-34;
