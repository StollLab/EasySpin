% gaussian  Gaussian line shape
%
%   ya = gaussian(x,x0,fwhm)
%   ya = gaussian(x,x0,fwhm,diff)
%   ya = gaussian(x,x0,fwhm,diff,phase)
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
%            2 second and so on, -1 the integral with -infinity
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

if nargin<3
  error('gaussian needs three inputs: x, x0, and fwhm.');
end
if nargin<4, diff = 0; end
if nargin<5, phase = 0; end

% Check arguments
%-------------------------------------------------------------------------------
if numel(x0)~=1 || ~isreal(x0)
  error('2nd input (x0) must be a real number.');
end

if numel(fwhm)~=1 || ~isreal(fwhm) || fwhm<=0
  error('3rd input (fwhm) must be a positive real number.');
end

if numel(diff)~=1 || ~isreal(diff) || mod(diff,1)~=0 || diff<-1
  error('4th input (diff) must be -1, 0, 1, 2, ...');
end

if numel(phase)~=1 || ~isreal(phase)
  error('5th input (phase) must be a real value.');
end

nonzeroPhase = mod(phase,2*pi)~=0;

if diff<0 && nonzeroPhase
  error('Cannot compute phased lineshape for integral.');
end

calcDispersion = nonzeroPhase || nargout>1;

% Compute Gaussian lineshape
%-------------------------------------------------------------------------------
% sig = distance from center to inflexion point (== standard deviation)
sig = fwhm/sqrt(2*log(2))/2;

if diff==-1
  
  yabs = 1/2*(1+erf((x-x0)/sig/sqrt(2)));
  ydisp = NaN(size(yabs));
  
else
  
  % absorption
  z = (x-x0)/sig;
  yabs = sqrt(2/pi)/2/sig*(-1/sig)^diff*2^(-diff/2)*...
    hermitepoly(z/sqrt(2),diff).*exp(-z.^2/2);
  
  % dispersion
  if calcDispersion
    k = z/sqrt(2);
    n = diff;
    ydisp = sqrt(2/pi)*2/sqrt(pi)/sig/2*(1/sqrt(2)/sig)^n * dawsonFderiv(k,n);
  end
  
end

% Phase rotation
if nonzeroPhase
  y1 =  yabs*cos(phase) + ydisp*sin(phase);
  y2 = -yabs*sin(phase) + ydisp*cos(phase);
else
  y1 = yabs;
  if calcDispersion
    y2 = ydisp;
  end
end

return
%===============================================================================

function y = dawsonFderiv(k,n)
y = (-1)^n * (hermitepoly(k,n).*dawsonF(k) - Gpoly(k,n-1));
return


function y = dawsonF(x)
y = real(sqrt(pi)/2i*(faddeeva(x,38)-exp(-x.^2)));
return

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
return

function y = hermitepoly(x,n)
% Calculates the Hermite polynomial of n-th order
% (probabilist definition, where highest order coefficient = 1 for all n)
switch n
  case 0
    y = 1;
  case 1
    y = 2*x;
  otherwise
    y = 2*x.*hermitepoly(x,n-1)-2*(n-1)*hermitepoly(x,n-2);
end
return
