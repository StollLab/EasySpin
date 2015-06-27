% eps0  Electric constant, vacuum permittivity
%
%   e = eps0
%   [e,sigma] = eps0
%
%   Returns the electric constant (also known as vacuum
%   permittivity) in SI units, farads per meter.
%
%   sigma is the standard uncertainty, which is zero for
%   this constant.

function [e,sigma] = eps0

e = 1/mu0/clight^2;
sigma = 0;
