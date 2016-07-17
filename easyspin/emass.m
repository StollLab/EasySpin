% emass  Electron mass 
%
%   me = emass
%   [me,sigma] = emass
%
%   Returns the mass of the electron in SI units, kg.
%   sigma is the standard uncertainty.

function [me,sigma] = emass

% 2014 CODATA value
% Concise form 	 9.109 383 56(11) x 10-31 kg 

me =    9.10938356e-31;
sigma = 0.00000011e-31;
