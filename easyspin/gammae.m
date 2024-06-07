% gammae  Electron gyromagnetic ratio
%
%   val = gammae
%   [val,sig] = gammae
%
%   Returns the value of the gyromagnetic ratio of the electron, in SI
%   units, rad s^-1 T^-1.
%   sigma is the standard deviation.

function [value,sigma] = gammae

% 2022 CODATA value, with added negative sign

value = -1.76085962784e11;
sigma =  0.00000000055e11;

end
