% hartree  Atomic unit of energy 
%
%   E = hartree
%   [E,sigma] = hartree
%
%   Returns the Hartree energy, the atomic unit of energy
%   in J. sigma is the standard uncertainty (2010 CODATA).

function [E,sigma] = hartree

% 2010 CODATA value
%Value 	 4.359 744 34 x 10-18 J
% Standard uncertainty 	 0.000 000 19 x 10-18 J
%  Relative standard uncertainty 	  4.4 x 10-8
% Concise form 	 4.359 744 34(19) x 10-18 J

E =     4.35974434e-18;
sigma = 0.00000019e-18;

return
