% mu0  Magnetic constant, vacuum permeability
%
%   m = mu0
%   [m,sigma] = mu0
%
%   Returns the magnetic constant (also known as vacuum
%   permeability) in SI units (N A^-2 = T^2 m^3 J^-1).
%
%   sigma is the standard uncertainty, which is zero for
%   this constant.

function [m,sigma] = mu0

m = 4*pi*1e-7;
sigma = 0;

return
