function ok = test(opt,olddata)

% Test against a few explicit values
%===============================================================================
% Explicit exressions from
%   M. A. Blanco, M. Flóres, M. Bermejo
%   Evaluation of the rotation matrices in the basis of real spherical harmonics
%   Journal of Molecular Structure (Theochem) 419 (1997) 19–27
%   https://doi.org/10.1016/S0166-1280(97)00185-1
%   Table 1
% See also
%   C.D.H. Chisholm
%   Group theoretical techniques in quantum chemistry
%   Academic Press, 1976

rng(3434);
n = 100;
theta = rand(1,n)*pi;
phi = rand(1,n)*2*pi;

LMf{1} = {0  0  @(t,p) ones(size(t))}; 
LMf{2} = {1 -1  @(t,p) sin(t).*sin(p)};
LMf{3} = {1  0  @(t,p) cos(t)};
LMf{4} = {1 +1  @(t,p) sin(t).*cos(p)};
LMf{5} = {2 -2  @(t,p) sqrt(3)/2*sin(t).^2.*sin(2*p)};
LMf{6} = {2 -1  @(t,p) sqrt(3)/2*sin(2*t).*sin(p)};
LMf{7} = {2  0  @(t,p) 1/2*(3*cos(t).^2-1)};
LMf{8} = {2 +1  @(t,p) sqrt(3)/2*sin(2*t).*cos(p)};
LMf{9} = {2 +2  @(t,p) sqrt(3)/2*sin(t).^2.*cos(2*p)};

for k = 1:numel(LMf)
  L = LMf{k}{1};
  M = LMf{k}{2};
  f = LMf{k}{3};
  vref = sqrt((2*L+1)/4/pi)*f(theta,phi);
  v = spherharm(L,M,theta,phi,'r');
  ok(k) = areequal(v,vref,1e-10,'rel');
end
