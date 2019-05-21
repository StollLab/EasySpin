% nmass   neutron mass 
%
%   mn = nmass
%   [mn,sigma] = nmass
%
%   Returns the mass of the neutron, in kg.
%   sigma is the standard uncertainty.

function [mn,sigma] = nmass

% 2018 CODATA value

mn =    1.67492749804e-27;
sigma = 0.00000000095e-27;

return
