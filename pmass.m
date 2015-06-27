% pmass  Proton mass 
%
%   mp = pmass
%   [mp,sigma] = pmass
%
%   Returns the mass of the proton, in kg.
%   sigma is the standard uncertainty.

function [mp,sigma] = pmass

% 2014 CODATA value
% Concise form 	 1.672 621 898(21) x 10-27 kg

mp =    1.672621898e-27;
sigma = 0.000000021e-27;
