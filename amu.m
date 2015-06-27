% amu  Atomic mass unit 
%
%   u = amu
%   [u,sigma] = amu
%
%   Returns the amotic mass unit
%   in kg. sigma is the standard
%   uncertainty.

function [u,sigma] = amu

% 2014 CODATA value

% Value 1.660 539 040 x 10-27 kg
% Standard uncertainty 0.000 000 020 x 10-27 kg
% Relative standard uncertainty 1.2 x 10-8
% Concise form 1.660 530 040(20) x 10-27 kg
u =     1.660539040e-27;
sigma = 0.000000020e-27;

return
