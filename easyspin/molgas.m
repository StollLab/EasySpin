% molgas  molar gas constant 
%
%   R = molgas
%
%   Returns the molar gas constant in SI units, joule per mole per kelvin.

function R = molgas

% 2018 CODATA value

R = avogadro*boltzm; % exact, as of SI redefinition 2019

return
