% nmagn  Nuclear magneton 
%
%   mu = nmagn
%   [mu,sigma] = nmagn
%
%   Returns the nuclear magneton
%   in Joule/Tesla. sigma is the
%   standard uncertainty (2006 CODATA).

function [mu_n,sigma] = nmagn

% 2006 CODATA recommended values

% Value 5.050 783 23 x 10-27 J T-1
% Standard uncertainty 0.000 000 13 x 10-27 J T-1
% Relative standard uncertainty 2.5 x 10-8
% Concise form 5.050 783 23(13) x 10-27 J T-1

mu_n =  5.05078323e-27;
sigma = 0.00000013e-27;

return
