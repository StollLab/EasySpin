% eps0  Electric constant, vacuum permittivity
%
%   e = eps0
%
%   Returns the electric constant (also known as vacuum
%   permittivity) in SI units, farads per meter.

function e = eps0

e = 1/mu0/clight^2;

return
