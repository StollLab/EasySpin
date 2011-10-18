% amu  Atomic mass unit 
%
%   u = amu
%   [u,sigma] = amu
%
%   Returns the amotic mass unit
%   in kg. sigma is the standard
%   uncertainty (2006 CODATA).

function [u,sigma] = amu

% 2006 CODATA value

% Value 1.660 538 782 x 10-27 kg
% Standard uncertainty 0.000 000 083 x 10-27 kg
% Relative standard uncertainty 5.0 x 10-8
% Concise form 1.660 538 782(83) x 10-27 kg
u =     1.660538782e-27;
sigma = 0.000000083e-27;

return
