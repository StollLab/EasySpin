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
%   - x:     abscissa vector
%   - x0:    center of the lineshape function
%   - fwhm:  full width at half height
%   - diff:  derivative. 0 is no derivative, 1 first,
%            2 second and so on, -1 the integral with -infinity
%            as lower limit. 0 is the default.
%   - phase: phase rotation, mixes absorption and dispersion.
%            phase=pi/2 puts dispersion signal into ya
%
%   Output:
%   - ya:   absorption function values for abscissa x
%   - yd:   dispersion function values for abscissa x

% For Gaussian dispersion equations, see
%   Pake and Purcell, Phys. Rev. 74(9), 1184-1188 (1948)
%     Line Shapes in Nuclear Paramagnetism
% For Derivatives of the Dawson function, see
%   R. Barakat, J. Quant. Spectrosc. Radiat. Transfer 11(11), 1729-1730, 1971.

function [yabs,ydisp] = gaussian(x,x0,fwhm,diff,phase)

if (nargin==0), help(mfilename); return; end

if (nargin<4), diff = 0; end
if (nargin<5), phase = 0; end

if any(fwhm<=0) || any(~isreal(fwhm))
  error('fwhm must be positive and real!');
end

if numel(fwhm)>1
  error('fwhm must contain 1 element!');
end

if any(diff<-1)
  error('Cannot compute lineshape for diff=%d.',diff);
end

DoPhase = 0;
if (nargin>4)
  DoPhase = 1;
  if numel(phase)>1
    error('phase must contain 1 element.');
  elseif ~isreal(phase)
    error('phase must be real.');
  elseif (diff<0)
    error('Cannot compute phased lineshape for integral.');
  end
end


% Compute Gaussian lineshape
%------------------------------------------------------------------
% gamma = distance from center to inflexion point
gamma = fwhm/sqrt(2*log(2))/2;

if (diff==-1)
  yabs = 1/2*(1+erf((x-x0)/gamma/sqrt(2)));
  ydisp = NaN;
else
  % absorption
  k = (x-x0)/gamma;
  yabs = sqrt(2/pi)/2/gamma*(-1/gamma)^diff*2^(-diff/2)*hermitepoly(k/sqrt(2),diff).*exp(-k.^2/2);

  % dispersion
  pre = sqrt(2/pi);
  if (nargout>1 || DoPhase) 
    k = (x-x0)/gamma/sqrt(2);
    n = diff;
    ydisp = pre*2/sqrt(pi)/gamma/2*(1/sqrt(2)/gamma)^n *...
      (-1)^n * (hermitepoly(k,n).*dawsonF(k) - G_Barakat(k,n-1));
  end
end

% Phase rotation
if (DoPhase)
  yabs1 =  yabs*cos(phase) + ydisp*sin(phase);
  ydisp = -yabs*sin(phase) + ydisp*cos(phase);
  yabs = yabs1;
end

return

function y = dawsonF(x)
y = real(sqrt(pi)/2i*(faddeeva(x,38)-exp(-x.^2)));
return

function y = G_Barakat(x,n)
% caclulates the G_n(x) as given in 
% R. Barakat, J. Quant. Spectrosc. Radiat. Transfer 11(11), 1729-1730, 1971.
switch n
  case -1
    y = 0;
  case 0
    y = 1;
  case 1
    y = 2*x;
  otherwise
    y = 2*x.*G_Barakat(x,n-1) - 2*n*G_Barakat(x,n-2);
end
return

function y = hermitepoly(x,n)
% calculates the Hermite polynomial of n-th order
% (probabilist definition, where highest order coefficient = 1 for all n)
switch n
  case 0, 
    y = 1;
  case 1
    y = 2*x;
  otherwise
    y = 2*x.*hermitepoly(x,n-1)-2*(n-1)*hermitepoly(x,n-2);
end
return



