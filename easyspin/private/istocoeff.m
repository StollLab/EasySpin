% istocoeff takes a 3x3 interaction tensor A as input and gives as output the
% coefficient for the isotropic ISTO component and the 5 coefficients
% corresponding to the anisotropic ISTO components
% Alternatively, it takes the 3 principal values.

function [F0,F1,F2] = istocoeffnew(A)

if numel(A)==1
  
  % Just isotropic value given
  
  F0 = -sqrt(1/3)*(3*A);
  F1 = zeros(3,1);
  F2 = zeros(5,1);
  
elseif numel(A)==3
  
  % Only 3 principal values given
  
  T0 = -sqrt(1/3)*[1 1 1];
  T2 = [1/2 -1/2 0; 0 0 0; sqrt(2/3)*[-1/2 -1/2 1]; 0 0 0; 1/2 -1/2 0];
  F0 = T0*A(:);
  F1 = zeros(3,1);
  F2 = T2*A(:);

elseif numel(A)==9
  
  % Full 3x3 tensor given
  
  x = 1;
  y = 2;
  z = 3;
  
  % Rank 0 (L = 0)
  %-------------------------------------
  F0 = -(1/sqrt(3))*(A(x,x)+A(y,y)+A(z,z)); % M = 0
  
  % Rank 1 (L = 1)
  %-------------------------------------
  F1 = zeros(3,1);
  F1(1) = -0.5*(A(z,x) - A(x,z) + 1i*(A(z,y) - A(y,z)));  % M = +1
  F1(3) = -0.5*(A(z,x) - A(x,z) - 1i*(A(z,y) - A(y,z)));  % M = -1
  F1(2) = -(1i/sqrt(2))*(A(x,y) - A(y,x));                % M =  0
  
  % Rank 2 (L = 2)
  %-------------------------------------
  F2 = zeros(5,1);
  F2(1) =  0.5*((A(x,x) - A(y,y)) + 1i*(A(x,y) + A(y,x)));  % M = +2
  F2(5) =  0.5*((A(x,x) - A(y,y)) - 1i*(A(x,y) + A(y,x)));  % M = -2
  F2(2) = -0.5*((A(x,z) + A(z,x)) + 1i*(A(y,z) + A(z,y)));  % M = +1
  F2(4) = +0.5*((A(x,z) + A(z,x)) - 1i*(A(y,z) + A(z,y)));  % M = -1
  F2(3) = sqrt(2/3) * (A(z,z)-0.5*(A(x,x) + A(y,y)));       % M =  0

else
  
  error('Cartesian tensor has wrong number of components.');

end

return