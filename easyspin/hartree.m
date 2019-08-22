% hartree  Atomic unit of energy 
%
%   E = hartree
%   [E,sigma] = hartree
%
%   Returns the Hartree energy, the atomic unit of energy
%   in joule. sigma is the standard uncertainty.

function [E,sigma] = hartree

% 2018 CODATA value

E =     4.3597447222071e-18;
sigma = 0.0000000000085e-18;

return
