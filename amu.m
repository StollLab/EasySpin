% amu  Atomic mass unit 
%
%   u = amu
%   [u,sigma] = amu
%
%   Returns the amotic mass unit
%   in kg. sigma is the standard
%   uncertainty (2010 CODATA).

function [u,sigma] = amu

% 2010 CODATA value

% Value 1.660 538 921 x 10-27 kg
% Standard uncertainty 0.000 000 073 x 10-27 kg
% Relative standard uncertainty 4.4 x 10-8
% Concise form 1.660 538 921(73) x 10-27 kg
u =     1.660538921e-27;
sigma = 0.000000073e-27;

return
