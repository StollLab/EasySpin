% molgas  molar gas constant 
%
%   R = molgas
%   [R,sigma] = molgas
%
%   Returns the molar gas constant in J/mol/K.
%   sigma is the standard uncertainty (2010 CODATA).

function [R,sigma] = molgas

% 2010 CODATA value

%8.314 4621 J mol-1 K-1
% Standard uncertainty 	 0.000 0075 J mol-1 K-1
%  Relative standard uncertainty 	  9.1 x 10-7
% Concise form 	 8.314 4621(75) J mol-1 K-1 
R =     8.3144621;
sigma = 0.0000075;

return
