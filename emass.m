% emass  Electron mass 
%
%   me = emass
%   [me,sigma] = emass
%
%   Returns the mass of the electron in kg.
%   sigma is the standard uncertainty (2010 CODATA).

function [me,sigma] = emass

% 2010 CODATA value
%Value 	 9.109 382 91 x 10-31 kg
% Standard uncertainty 	 0.000 000 40 x 10-31 kg
%  Relative standard uncertainty 	  4.4 x 10-8
% Concise form 	 9.109 382 91(40) x 10-31 kg 

me =    9.10938291e-31;
sigma = 0.00000040e-31;
