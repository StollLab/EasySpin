% pmass  Proton mass 
%
%   mp = pmass
%   [mp,sigma] = pmass
%
%   Returns the mass of the proton
%   in kg. sigma is the standard
%   uncertainty (2010 CODATA).

function [mp,sigma] = pmass

% 2010 CODATA value

% Value 	 1.672 621 777 x 10-27 kg
% Standard uncertainty 	 0.000 000 074 x 10-27 kg
%  Relative standard uncertainty 	  4.4 x 10-8
% Concise form 	 1.672 621 777(74) x 10-27 kg

mp =    1.672621777e-27;
sigma = 0.000000074e-27;

return
