% istocoeff  Irreducbile spherical tensors from cartesian tensor
%
% istocoeff converts a cartesian 3x3 interaction tensor matrix into
% its three irreducible spherical tensors (rank 0, 1, and 2).
%
% Input:
%    Fc  ... 3x3 cartesian interaction matrix (real-valued)
%            or 3x1 array of principal values
%            or single number if isotropic
%
% Output:
%    F0  ... rank-0 irreducbile spherical tensor (scalar)
%    F1  ... rank-1 irreducible spherical tensor (3x1 array)
%    F2  ... rank-2 irreducible spherical tensor (5x1 array)
%
% Reference:
%   Michael Mehring
%   Principles of High Resolution NMR in Solids, 2nd edition
%   Wiley, 1983
%   Appendix A, starting on p.288

function [F0,F1,F2] = istocoeff(Fc)

if ~isnumeric(Fc) || ~isreal(Fc)
  error('Input must be a real-valued scalar, a 3x1 array or a 3x3 array.');
end

% Initialize with all-zero
F0 = 0;   %#ok
F1 = zeros(3,1);
F2 = zeros(5,1);

% Array indices
x = 1;
y = 2;
z = 3;

if numel(Fc)==1
  
  % Just isotropic value given
  
  F0 = -sqrt(1/3)*(3*Fc);    % (0,0)
  
elseif numel(Fc)==3
  
  % Three principal values given  
  
  F0    = -sqrt(1/3)*(Fc(x)+Fc(y)+Fc(z));           % (0,0)
  F2(1) = +0.5*(Fc(x) - Fc(y));                     % (2,+2)
  F2(5) = +0.5*(Fc(x) - Fc(y));                     % (2,-2)
  F2(3) = sqrt(2/3) * (Fc(z)-0.5*(Fc(x) + Fc(y)));  % (2,0)
  
elseif numel(Fc)==9
  
  % Full 3x3 tensor given
  
  F0    = -sqrt(1/3)*(Fc(x,x)+Fc(y,y)+Fc(z,z ));                % (0,0)
  F1(1) = -0.5*(Fc(z,x) - Fc(x,z) + 1i*(Fc(z,y) - Fc(y,z)));    % (1,+1)
  F1(3) = -0.5*(Fc(z,x) - Fc(x,z) - 1i*(Fc(z,y) - Fc(y,z)));    % (1,-1)
  F1(2) = -(1i/sqrt(2))*(Fc(x,y) - Fc(y,x));                    % (1,0)
  F2(1) =  0.5*((Fc(x,x) - Fc(y,y)) + 1i*(Fc(x,y) + Fc(y,x)));  % (2,+2)
  F2(5) =  0.5*((Fc(x,x) - Fc(y,y)) - 1i*(Fc(x,y) + Fc(y,x)));  % (2,-2)
  F2(2) = -0.5*((Fc(x,z) + Fc(z,x)) + 1i*(Fc(y,z) + Fc(z,y)));  % (2,+1)
  F2(4) = +0.5*((Fc(x,z) + Fc(z,x)) - 1i*(Fc(y,z) + Fc(z,y)));  % (2,-1)
  F2(3) = sqrt(2/3) * (Fc(z,z)-0.5*(Fc(x,x) + Fc(y,y)));        % (2,0)

else
  
  error('Cartesian tensor has wrong size. It needs to be 1x1, 3x1, 1x3, or 3x3.');

end

return
