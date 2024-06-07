% faraday  Faraday constant 
%
%   F = faraday
%
%   Returns the Faraday constant in SI units, C mol^-1.

function F = faraday

% 2022 CODATA value
  
F = evolt*avogadro; % exact, as of 2019 SI redefinition

end
