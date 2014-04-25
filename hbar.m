% hbar    Planck constant over 2 pi
%
%   value = hbar
%   [value,sigma] = hbar
%
%   Returns the Planck constant over 2 pi, h/(2*pi), in Joule*second.
%   sigma is the standard uncertainty (2010 CODATA).

function [value,sigma] = hbar

% 2010 CODATA recommended value

% Value 	 1.054 571 726 x 10-34 J s
% Standard uncertainty 	 0.000 000 047 x 10-34 J s
%  Relative standard uncertainty 	  4.4 x 10-8
% Concise form 	 1.054 571 726(47) x 10-34 J s    
 
value = 1.054571726e-34;
sigma = 0.000000047e-34;

return
