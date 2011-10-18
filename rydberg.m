% rydberg  Rydberg constant 
%
%   Rinf = rydberg
%   [Rinf,sigma] = rydberg
%
%   Returns the Rydberg constant
%   in m^-1. sigma is the standard
%   uncertainty (2006 CODATA).

function [Rinf,sigma] = rydberg

% 2006 CODATA value

% Value 10 973 731.568 527 m-1
% Standard uncertainty 0.000 073 m-1
% Relative standard uncertainty 6.6 x 10-12
% Concise form 10 973 731.568 527(73) m-1

Rinf = 10973731.568527;
sigma =       0.000073;

return
