% faraday  Faraday constant 
%
%   F = faraday
%   [F,sigma] = faraday
%
%   Returns the Faraday constant in
%   C mol^-1. sigma is the standard
%   uncertainty (2006 CODATA).

function [F,sigma] = faraday

% 2006 CODATA value

% Value 96 485.3399 C mol-1
% Standard uncertainty 0.0024 C mol-1
% Relative standard uncertainty 2.5 x 10-8
% Concise form 96 485.3399(24) C mol-1

F = 96485.3399;
sigma = 0.0024;

return
