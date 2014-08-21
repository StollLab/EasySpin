% nmagn  Nuclear magneton 
%
%   mu = nmagn
%   [mu,sigma] = nmagn
%
%   Returns the nuclear magneton in SI units,
%   joule per tesla. sigma is the standard
%   uncertainty (2010 CODATA).

function [mu_n,sigma] = nmagn

% 2010 CODATA recommended values

% Value 	 5.050 783 53 x 10-27 J T-1
% Standard uncertainty 	 0.000 000 11 x 10-27 J T-1
%  Relative standard uncertainty 	  2.2 x 10-8
% Concise form 	 5.050 783 53(11) x 10-27 J T-1
mu_n =  5.05078353e-27;
sigma = 0.00000011e-27;

return
