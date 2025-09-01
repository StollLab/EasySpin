% gamman  Nuclear gyromagnetic ratio
%
%   val = gammae(Iso)
%
%   Returns the value of the gyromagnetic ratio for the nuclear isotopes listed
%   in Iso, in SI units (rad s^-1 T^-1).
%
%   Example:
%     g1H = gamman('1H')

function value = gamman(Iso)

value = nmagn*nucgval(Iso)/hbar;

end
