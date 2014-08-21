% planck  Planck constant 
%
%   h = planck
%   [h,sigma] = planck
%
%   Returns the Planck constant in SI units, joule times second.
%   sigma is the standard uncertainty (2010 CODATA).

function [h,sigma] = planck

% 2010 CODATA recommended value

%Value 	 6.626 069 57 x 10-34 J s
% Standard uncertainty 	 0.000 000 29 x 10-34 J s
%  Relative standard uncertainty 	  4.4 x 10-8
% Concise form 	 6.626 069 57(29) x 10-34 J s

h =     6.62606957e-34;
sigma = 0.00000029e-34;

return
