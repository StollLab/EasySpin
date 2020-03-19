% tensor_cart2sph  Irreducbile spherical tensors from rank-2 cartesian tensor
%
%    [T0,T1,T2] = tensor_cart2sph(Tc)
%
% tensor_cart2sph converts a cartesian tensor (3x3 matrix) into
% its three irreducible spherical tensors (ranks 0, 1, and 2).
%
% Input:
%    Tc  ... cartesian tensor (real-valued); 3x3 matrix, or 3x1 array of
%            principal values, or single number if isotropic
%
% Output:
%    T0  ... rank-0 irreducbile spherical tensor (scalar)
%    T1  ... rank-1 irreducible spherical tensor (3x1 array)
%    T2  ... rank-2 irreducible spherical tensor (5x1 array)

% Reference:
%   Michael Mehring
%   Principles of High Resolution NMR in Solids, 2nd edition
%   Wiley, 1983
%   Appendix A, starting on p.288

function [T0,T1,T2] = tensor_cart2sph(Tc)

if ~isnumeric(Tc) || ~isreal(Tc)
  error('Input must be a real-valued scalar, a 3x1 array or a 3x3 array.');
end

% Initialize with all-zero
T0 = 0;   %#ok
T1 = zeros(3,1);
T2 = zeros(5,1);

% Array indices
x = 1;
y = 2;
z = 3;

if numel(Tc)==1
  
  % Just isotropic value given
  
  T0 = -sqrt(1/3)*(3*Tc);    % (0,0)
  
elseif numel(Tc)==3
  
  % Three principal values given  
  
  T0    = -sqrt(1/3)*(Tc(x)+Tc(y)+Tc(z));           % (0,0)
  T2(1) = +0.5*(Tc(x) - Tc(y));                     % (2,+2)
  T2(5) = +0.5*(Tc(x) - Tc(y));                     % (2,-2)
  T2(3) = sqrt(2/3) * (Tc(z)-0.5*(Tc(x) + Tc(y)));  % (2,0)
  
elseif numel(Tc)==9
  
  % Full 3x3 tensor given
  
  T0    = -sqrt(1/3)*(Tc(x,x)+Tc(y,y)+Tc(z,z ));                % (0,0)
  T1(1) = -0.5*(Tc(z,x) - Tc(x,z) + 1i*(Tc(z,y) - Tc(y,z)));    % (1,+1)
  T1(3) = -0.5*(Tc(z,x) - Tc(x,z) - 1i*(Tc(z,y) - Tc(y,z)));    % (1,-1)
  T1(2) = -(1i/sqrt(2))*(Tc(x,y) - Tc(y,x));                    % (1,0)
  T2(1) =  0.5*((Tc(x,x) - Tc(y,y)) + 1i*(Tc(x,y) + Tc(y,x)));  % (2,+2)
  T2(5) =  0.5*((Tc(x,x) - Tc(y,y)) - 1i*(Tc(x,y) + Tc(y,x)));  % (2,-2)
  T2(2) = -0.5*((Tc(x,z) + Tc(z,x)) + 1i*(Tc(y,z) + Tc(z,y)));  % (2,+1)
  T2(4) = +0.5*((Tc(x,z) + Tc(z,x)) - 1i*(Tc(y,z) + Tc(z,y)));  % (2,-1)
  T2(3) = sqrt(2/3) * (Tc(z,z)-0.5*(Tc(x,x) + Tc(y,y)));        % (2,0)

else
  
  error('Cartesian tensor has wrong size. It needs to be 1x1, 3x1, 1x3, or 3x3.');

end

return
