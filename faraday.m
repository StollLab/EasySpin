% faraday  Faraday constant 
%
%   F = faraday
%   [F,sigma] = faraday
%
%   Returns the Faraday constant in
%   C mol^-1. sigma is the standard
%   uncertainty (2010 CODATA).

function [F,sigma] = faraday

% 2010 CODATA value
%Value 	 96 485.3365 C mol-1
% Standard uncertainty 	      0.0021 C mol-1
%  Relative standard uncertainty 	  2.2 x 10-8
% Concise form 	 96 485.3365(21) C mol-1 
 
F = 96485.3365;
sigma = 0.0021;

return
