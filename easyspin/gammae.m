% gammae  Electron gyromagnetic ratio
%
%   val = gammae
%   [val,sig] = gammae
%
%   Returns the value of the gyromagnetic ratio of the electron, in SI
%   units (radians per second per tesla, rad s^-1 T^-1).
%   sigma is the standard deviation.

function [value,sigma] = gammae

% CODATA 2018 value, with added negative sign

value = -1.76085963023e11;
sigma =  0.00000000053e11;

end
