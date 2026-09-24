% gaussian  Gaussian line shape
%
%   ya = gaussian(x,x0,fwhm)
%   ya = gaussian(x,x0,fwhm,diff)
%   y = gaussian(x,x0,fwhm,diff,phase)
%   [ya,yd] = gaussian(...)
%
%   Computes area-normalized Gaussian absorption and dispersion
%   line shapes or their derivatives.
%
%   Input:
%   - x:     abscissa vector (any units)
%   - x0:    center of the lineshape function (same units as x)
%   - fwhm:  full width at half height (same units as x)
%   - diff:  derivative. 0 is no derivative, 1 first,
%            2 second derivative, etc; -1 the integral with -infinity
%            as lower limit. 0 is the default.
%   - phase: phase rotation (radians), mixes absorption and dispersion.
%            phase=pi/2 puts dispersion signal into ya
%
%   Output:
%   - ya:   absorption function values for abscissa x
%   - yd:   dispersion function values for abscissa x

% For Gaussian dispersion equations, see
%   G. E. Pake and E. M. Purcell
%   Line Shapes in Nuclear Paramagnetism
%   Phys. Rev. 74(9), 1184-1188 (1948)
%   https://doi.org/10.1103/PhysRev.74.1184
% For derivatives of the Dawson function, see
%   R. Barakat
%   The derivatives of Dawson's function
%   J. Quant. Spectrosc. Radiat. Transfer 11(11), 1729-1730, 1971.
%   https://doi.org/10.1016/0022-4073(71)90151-8

function [y1,y2] = gaussian(x,x0,fwhm,diff,phase)

if nargin==0, help(mfilename); return; end

% Check number of input arguments
%-------------------------------------------------------------------------------
if nargin<3
  error('gaussian() needs at least three inputs: x, x0, and fwhm.');
end

if nargin<4, diff = 0; end
if nargin<5, phase = 0; end

% Check input arguments
%-------------------------------------------------------------------------------
if ~isscalar(x0) || ~isreal(x0)
  error('2nd input (x0, center) must be a real number.');
end

if ~isscalar(fwhm) || ~isreal(fwhm) || fwhm<=0
  error('3rd input (fwhm, full width at half maximum) must be a real and positive number.');
end

if ~isscalar(diff) || diff<-1 || mod(diff,1)
  error('4th input (diff, lineshape derivative) must be -1, 0, 1, 2, etc.');
end

if ~isscalar(phase) || ~isreal(phase)
  error('5th input (phase) must be a real number.');
end

nonzeroPhase = mod(phase,2*pi)~=0;

returnDispersion = nargout>1;
calcDispersion = nonzeroPhase || returnDispersion;

% Compute Gaussian lineshape
%-------------------------------------------------------------------------------
% sig = distance from center to inflexion point (== standard deviation)
sig = fwhm/sqrt(2*log(2))/2;
k = (x-x0)/sig/sqrt(2);

if diff==-1
  
  % absorption
  yabs = 1/2*(1+erf(k));

  % dispersion
  if calcDispersion
    ydisp = sqrt(2/pi)*2/sqrt(pi)/sig/2 * ...  % prefactor
      sig*sqrt(2) *...  % dx/dk
      dawsonFintegral(k);  % integral of Dawson function
  end
  
else
  
  % absorption
  yabs = sqrt(2/pi)/2/sig*(-1/sig)^diff*2^(-diff/2) * ...
    hermitepoly(k,diff).*exp(-k.^2);
  
  % dispersion
  if calcDispersion
    ydisp = sqrt(2/pi)*2/sqrt(pi)/sig/2*(1/sqrt(2)/sig)^diff * ...
      dawsonFderiv(k,diff);
  end
  
end

% Phase rotation
%-------------------------------------------------------------------------------
if nonzeroPhase
  y1 =  yabs*cos(phase) + ydisp*sin(phase);
  y2 = -yabs*sin(phase) + ydisp*cos(phase);
else
  y1 = yabs;
  if calcDispersion
    y2 = ydisp;
  end
end

end
%===============================================================================

function y = dawsonF(x)
% Calculates the Dawson function F(x) (one-sided Fourier-Laplace sine transform
% of the Gaussian function)
y = sqrt(pi)/2i * (faddeeva(x)-exp(-x.^2));
y = real(y);
end

function y = dawsonFderiv(x,n)
% Calculates the n-th derivative of the Dawson function, d^F(x)/dx^n
y = (-1)^n * (hermitepoly(x,n).*dawsonF(x) - Gpoly(x,n-1));
end

function y = dawsonFintegral(x)
% Calculates the integral of the Dawson function from 0 to x,
%   x^2/2*2F2([1,1],[3/2,2],-x^2)
% (even function of x)
y = zeros(size(x));
ax = abs(x);

% Small arguments: Gauss-Legendre quadrature over [0,|x|]
xmax = 6;
small = ax<=xmax;
if any(small(:))
  % Golub-Welsch algorithm for nodes u and weights w over [-1,1]
  n = 40;
  beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
  [V,D] = eig(diag(beta,1)+diag(beta,-1));
  [u,idx] = sort(diag(D));
  w = 2*V(1,idx).^2;
  a = ax(small);
  t = a(:)/2*(u.'+1);
  y(small) = a(:)/2.*(dawsonF(t)*w.');
end

% Large arguments: asymptotic expansion
%   (1/2)*log(2x) + gamma/4 - sum_k (2k-1)!!/2^k/(4k)/x^(2k)
if any(~small(:))
  a = ax(~small);
  EulerGamma = 0.57721566490153286;
  y_ = log(2*a)/2 + EulerGamma/4;
  c = 1;
  for k = 1:30
    c = c*(2*k-1)/2;  % (2k-1)!!/2^k
    y_ = y_ - c/(4*k)./a.^(2*k);
  end
  y(~small) = y_;
end
end

function y = Gpoly(x,n)
% Calculates the G_n(x) polynomial as given in
%   R. Barakat
%   The derivatives of Dawson's function
%   J. Quant. Spectrosc. Radiat. Transfer 11(11), 1729-1730, 1971.
%   https://doi.org/10.1016/0022-4073(71)90151-8
switch n
  case -1
    y = 0;
  case 0
    y = 1;
  case 1
    y = 2*x;
  otherwise
    y = 2*x.*Gpoly(x,n-1) - 2*n*Gpoly(x,n-2);
end
end

function y = hermitepoly(x,n)
% Calculates the Hermite polynomial of n-th order
% (physics definition, where highest-order coefficient = 2^n)
switch n
  case 0
    y = 1;
  case 1
    y = 2*x;
  otherwise
    y = 2*x.*hermitepoly(x,n-1) - 2*(n-1)*hermitepoly(x,n-2);
end
end
