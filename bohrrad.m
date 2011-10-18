% bohrrad  Bohr radius 
%
%   r = bohrrad
%   [r,sigma] = bohrrad
%
%   Returns the Bohr radius in units of meters.
%   sigma is the standard uncertainty (2006 CODATA).

function [r,sigma] = bohrrad

% 2006 CODATA value

% Value 0.529 177 208 59 x 10-10 m
% Standard uncertainty 0.000 000 000 36 x 10-10 m
% Relative standard uncertainty 6.8 x 10-10
% Concise form 0.529 177 208 59(36) x 10-10 m

r =     0.52917720859e-10;
sigma = 0.00000000036e-10;

return
