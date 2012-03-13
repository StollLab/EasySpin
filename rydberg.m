% rydberg  Rydberg constant 
%
%   Rinf = rydberg
%   [Rinf,sigma] = rydberg
%
%   Returns the Rydberg constant
%   in m^-1. sigma is the standard
%   uncertainty (2010 CODATA).

function [Rinf,sigma] = rydberg

% 2010 CODATA value

% Value 	 10 973 731.568 539 m-1
% Standard uncertainty 	          0.000 055 m-1
%  Relative standard uncertainty 	  5.0 x 10-12
% Concise form 	 10 973 731.568 539(55) m-1 

Rinf = 10973731.568539;
sigma =       0.000055;

return
