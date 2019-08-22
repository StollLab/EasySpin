% faraday  Faraday constant 
%
%   F = faraday
%
%   Returns the Faraday constant in SI units, coulomb per mole.

function F = faraday

% 2018 CODATA value
  
F = evolt*avogadro; % exact, as of 2019 SI redefinition

return
