% avogadro  Avogadro number 
%
%   N = avogadro
%   [N,sigma] = avogadro
%
%   Returns the Avogadro constant
%   in mol^-1. sigma is the standard
%   uncertainty.

function [N_A,sigma] = avogadro

% 2014 CODATA value
% Value 6.022 140 857 x 10^23 mol-1
% Standard uncertainty 0.000 000 074 x 10^23 mol-1
% Relative standard uncertainty 1.2 x 10^-8
% Concise form 6.022 140 857(74) x 10^23 mol-1 

N_A =   6.022140857e23;
sigma = 0.000000074e23;

return
