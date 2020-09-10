% spherharm    Complex- and real-valued spherical harmonics
%
%   y = spherharm(L,M,theta,phi)
%   y = spherharm(L,M,theta,phi,'r')
%
%   Computes spherical harmonic Y_L,M(theta,phi) with L>=0 and -L<=M<=L.
%
%   theta is the angle down from the z axis (colatitude), and phi is the
%   counterclockwise angle off the x axis in the xy plane (longitude).
%   theta and phi can be scalars or arrays (of the same size). Both angles
%   are assumed to be in units of radians.
%
%   If the option 'r' is included, the real-valued spherical harmonics
%   are returned, using cos(M*phi) for M>=0, and sin(abs(M)*phi) for M<0.
%
%   The complex-valued spherical harmonics evaluated by spherharm() include the
%   Condon-Shortley phase (-1)^M.
%
%   The sign of the real-valued harmonics is defined such that they give
%   nonnegative values near theta=0 and phi=0 for all L and M.

% The real-valued spherical harmonics are defined according to
%   M. A. Blanco, M. Flóres, M. Bermejo
%   Evaluation of the rotation matrices in the basis of real spherical harmonics
%   Journal of Molecular Structure (Theochem) 419 (1997) 19–27
%   https://doi.org/10.1016/S0166-1280(97)00185-1
%   Table 1
% See also
%   C.D.H. Chisholm
%   Group theoretical techniques in quantum chemistry
%   Academic Press, 1976

function y = spherharm(L,M,theta,phi,Type)

if nargin==0, help(mfilename); return; end

if nargin<4
  error('At least four inputs are required: L, M, theta, phi.');
end

if nargin>5
  error('At most five inputs are allowed: L, M, theta, phi, ''r''.');
end

if nargin<5, Type = 'c'; end

if numel(L)~=1 || L<0 || mod(L,1)
  error('L (first input) must be a non-negative integer (0, 1, 2, etc).');
end

if numel(M)~=1 || abs(M)>L || mod(M,1)
  error('M (second input) must be an integer between -L and L.');
end

if numel(theta)~=numel(phi)
  error('theta and phi (3rd and 4th input) must have the same size.');
end
if ~isreal(theta)
  error('theta (3rd input) must be real-valued.')
end
if ~isreal(phi)
  error('phi (4th input) must be real-valued.')
end

if numel(Type)~=1 || ~ischar(Type) || ~any(Type=='rc')
  error('If given, fourth input must be ''r'' or ''c''.');
end
realHarmonics = Type=='r';

if realHarmonics
  % do not use CS phase, to make all explicit expressions unsigned
  N = normfactor(L,abs(M));
  P = plegendre(L,abs(M),cos(theta),false);
  if M>0
    y = sqrt(2) * N * P.* cos(M*phi);
  elseif M<0
    y = sqrt(2) * N * P.* sin(abs(M)*phi);
  else
    y = N * P;
  end
else
  y = normfactor(L,M) * plegendre(L,M,cos(theta),true) .* exp(1i*M*phi);
end

end
%===============================================================================

function N = normfactor(L,M)
% Calculate normalization factor
  p = prod(L-abs(M)+1:L+abs(M));
  if M>0, p = 1/p; end
  N = sqrt((2*L+1)/(4*pi)*p);
end
