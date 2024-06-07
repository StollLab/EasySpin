% molgas  molar gas constant 
%
%   R = molgas
%
%   Returns the molar gas constant in SI units, J mol^-1 K^-1

function R = molgas

% 2022 CODATA value

R = avogadro*boltzm; % exact, as of SI redefinition 2019

end
