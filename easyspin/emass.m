% emass  Electron mass 
%
%   me = emass
%   [me,sigma] = emass
%
%   Returns the mass of the electron in SI units, kg.
%   sigma is the standard uncertainty.

function [me,sigma] = emass

% 2018 CODATA value

me =    9.1093837015e-31;
sigma = 0.0000000028e-31;

return
