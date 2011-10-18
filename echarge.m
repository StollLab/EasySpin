% echarge  Elementary charge 
%
%   e = echarge
%   [e,sigma] = echarge
%
%   Returns the elementary charge
%   in C. sigma is the standard
%   uncertainty (2006 CODATA).

function [e,sigma] = echarge

% 2006 CODATA value
% Value 1.602 176 487 x 10-19 C
% Standard uncertainty 0.000 000 040 x 10-19 C
% Relative standard uncertainty 2.5 x 10-8
% Concise form 1.602 176 487(40) x 10-19 C

e =     1.602176487e-19;
sigma = 0.000000040e-19;

return
