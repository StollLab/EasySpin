% pmass  Proton mass 
%
%   mp = pmass
%   [mp,sigma] = pmass
%
%   Returns the mass of the proton
%   in kg. sigma is the standard
%   uncertainty (2006 CODATA).

function [mp,sigma] = pmass

% 2006 CODATA value

% Value 1.672 621 637 x 10-27 kg
% Standard uncertainty 0.000 000 083 x 10-27 kg
% Relative standard uncertainty 5.0 x 10-8
% Concise form 1.672 621 637(83) x 10-27 kg

mp =    1.672621637e-27;
sigma = 0.000000083e-27;

return
