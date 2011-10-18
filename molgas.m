% molgas  molar gas constant 
%
%   R = molgas
%   [R,sigma] = molgas
%
%   Returns the molar gas constant in J/mol/K.
%   sigma is the standard uncertainty (2006 CODATA).

function [R,sigma] = molgas

% 2006 CODATA value

% Value 8.314 472 J mol-1 K-1
% Standard uncertainty 0.000 015 J mol-1 K-1
% Relative standard uncertainty 1.7 x 10-6
% Concise form 8.314 472(15) J mol-1 K-1

R =     8.314472;
sigma = 0.000015;

return
