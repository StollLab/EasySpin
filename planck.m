% planck  Planck constant 
%
%   h = planck
%   [h,sigma] = planck
%
%   Returns the Planck constant in Joule*second.
%   sigma is the standard uncertainty (2006 CODATA).

function [h,sigma] = planck

% 2006 CODATA recommended value

% Value 6.626 068 96 x 10-34 J s
% Standard uncertainty 0.000 000 33 x 10-34 J s
% Relative standard uncertainty 5.0 x 10-8
% Concise form 6.626 068 96(33) x 10-34 J s

h =     6.62606896e-34;
sigma = 0.00000033e-34;

return
