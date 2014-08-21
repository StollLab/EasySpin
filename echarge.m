% echarge  Elementary charge 
%
%   e = echarge
%   [e,sigma] = echarge
%
%   Returns the elementary charge in SI units,
%   coulomb. sigma is the standard
%   uncertainty (2010 CODATA).

function [e,sigma] = echarge

% 2010 CODATA value
%Value 	 1.602 176 565 x 10-19 C
% Standard uncertainty 	 0.000 000 035 x 10-19 C
%  Relative standard uncertainty 	  2.2 x 10-8
% Concise form 	 1.602 176 565(35) x 10-19 C
 
e =     1.602176565e-19;
sigma = 0.000000035e-19;

return
