% nmass   neutron mass 
%
%   mn = nmass
%   [mn,sigma] = nmass
%
%   Returns the mass of the neutron, in kg.
%   sigma is the standard uncertainty.

function [mn,sigma] = nmass

% 2014 CODATA value
% Concise form 	 1.674 927 471(21) x 10-27 kg

mn =    1.674927471e-27;
sigma = 0.000000021e-27;
