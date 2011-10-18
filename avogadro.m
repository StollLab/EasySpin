% avogadro  Avogadro number 
%
%   N = avogadro
%   [N,sigma] = avogadro
%
%   Returns the Avogadro constant
%   in mol^-1. sigma is the standard
%   uncertainty (2006 CODATA).

function [N_A,sigma] = avogadro

% 2006 CODATA value
% Value 6.022 141 579 x 1023 mol-1
% Standard uncertainty 0.000 000 30 x 1023 mol-1
% Relative standard uncertainty 5.0 x 10-8
% Concise form 6.022 141 79(30) x 1023 mol-1

N_A =   6.02214179e23;
sigma = 0.00000030e23;

return
