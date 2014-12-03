% istocoeff takes an interaction tensor as input and gives as output the
% coefficient for the isotropic ISTO component and the 5 coefficients
% corresponding to the anisotropic ISTO components

function [F0,F1,F2] = istocoeff(A)

x = 1;
y = 2;
z = 3;


% L = 0
%-------------------------------------
F0 = -(1/sqrt(3))*(A(x,x)+A(y,y)+A(z,z)); % M = 0

% L = 1
%-------------------------------------
F1 = zeros(1,3);
F1(1) = -0.5*(A(z,x) - A(x,z) + 1i*(A(z,y) - A(y,z)));  % M = +1
F1(3) = -0.5*(A(z,x) - A(x,z) - 1i*(A(z,y) - A(y,z)));  % M = -1
F1(2) = -(1i/sqrt(2))*(A(x,y) - A(y,x));                % M =  0

% L = 2
%-------------------------------------
F2 = zeros(1,5);
F2(1) =  0.5*((A(x,x) - A(y,y)) + 1i*(A(x,y) + A(y,x)));  % M = +2
F2(5) =  0.5*((A(x,x) - A(y,y)) - 1i*(A(x,y) + A(y,x)));  % M = -2
F2(2) = -0.5*((A(x,z) + A(z,x)) + 1i*(A(y,z) + A(z,y)));  % M = +1
F2(4) = +0.5*((A(x,z) + A(z,x)) - 1i*(A(y,z) + A(z,y)));  % M = -1
F2(3) = sqrt(2/3) * (A(z,z)-0.5*(A(x,x) + A(y,y)));       % M =  0

return