% emass  Electron mass 
%
%   me = emass
%   [me,sigma] = emass
%
%   Returns the mass of the electron in kg.
%   sigma is the standard uncertainty (2006 CODATA).

function [me,sigma] = emass

% 2006 CODATA value
% Value 9.109 382 15 x 10-31 kg
% Standard uncertainty 0.000 000 45 x 10-31 kg
% Relative standard uncertainty 5.0 x 10-8
% Concise form 9.109 382 15(45) x 10-31 kg

me =    9.10938215e-31;
sigma = 0.00000045e-31;

