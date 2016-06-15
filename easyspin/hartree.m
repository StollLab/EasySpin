% hartree  Atomic unit of energy 
%
%   E = hartree
%   [E,sigma] = hartree
%
%   Returns the Hartree energy, the atomic unit of energy
%   in joule. sigma is the standard uncertainty.

function [E,sigma] = hartree

% 2014 CODATA value
% Concise form 	 4.359 744 650(54) x 10-18 J

E =     4.359744650e-18;
sigma = 0.000000054e-18;

return
