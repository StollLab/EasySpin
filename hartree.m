% hartree  Atomic unif of energy 
%
%   E = hartree
%   [E,sigma] = hartree
%
%   Returns the atomic unit of energy
%   in J. sigma is the standard
%   uncertainty (2006 CODATA).

function [E,sigma] = hartree

% 2006 CODATA value
% Value 4.359 743 94 x 10-18 J
% Standard uncertainty 0.000 000 22 x 10-18 J
% Relative standard uncertainty 5.0 x 10-8
% Concise form 4.359 743 94(22) x 10-18 J

E = 4.35974394e-18;
sigma = 0.00000022e-18;

return
