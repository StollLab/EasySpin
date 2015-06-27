% bohrrad  Bohr radius 
%
%   r = bohrrad
%   [r,sigma] = bohrrad
%
%   Returns the Bohr radius in units of meters.
%   sigma is the standard uncertainty.

function [r,sigma] = bohrrad

% 2014 CODATA value

% Value 	 0.529 177 210 92 x 10-10 m
% Standard uncertainty 	 0.000 000 000 17 x 10-10 m
% Relative standard uncertainty 	  3.2 x 10-10
% Concise form 	 0.529 177 210 67(12) x 10-10 m
 
r =     0.52917721067e-10;
sigma = 0.00000000012e-10;

return
