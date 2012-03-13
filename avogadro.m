% avogadro  Avogadro number 
%
%   N = avogadro
%   [N,sigma] = avogadro
%
%   Returns the Avogadro constant
%   in mol^-1. sigma is the standard
%   uncertainty (2010 CODATA).

function [N_A,sigma] = avogadro

% 2010 CODATA value 6.022 141 29(27) x 10^23 mol-1

N_A =   6.02214129e23;
sigma = 0.00000027e23;

return
