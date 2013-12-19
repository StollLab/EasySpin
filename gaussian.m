% gaussian  Gaussian line shape 
%
%   ya = gaussian(x,x0,fwhm)
%   ya = gaussian(x,x0,fwhm,diff)
%   [ya,yd] = ...
%
%   Returns a area-normalized Gaussian line shape
%   or one of its derivatives.
%
%   Input:
%   - x:    Abscissa vector
%   - x0:   Centre of the lineshape function
%   - fwhm: Full width at half height
%   - diff: Derivative. 0 is no derivative, 1 first,
%           2 second, -1 the integral with -infinity
%           as lower limit. 0 is the default.
%
%   Output:
%   - ya:   absorption function values for abscissa x
%   - yd:   dispersion function values for abscissa x

% For Gaussian dispersion equations, see
% Pake and Purcell, Phys. Rev. 74(9), 1184-1188 (1948)
% Line Shapes in Nuclear Paramagnetism

function [yabs,ydisp] = gaussian(x,x0,fwhm,diff)

if (nargin==0), help(mfilename); return; end

if (nargin==3), diff = 0; end

if any(fwhm<=0) || any(~isreal(fwhm))
  error('fwhm must be positive and real!');
end

if numel(fwhm)>1
  error('fwhm must contain 1 element!');
end

if any(diff<-1) || any(diff>2)
  error('Cannot compute lineshape for derivative %d.',diff);
end

% Compute Gaussian lineshape
%------------------------------------------------------------------
% gamma = distance from center to inflexion point
gamma = fwhm/sqrt(2*log(2));
pre = sqrt(2/pi);
k = (x-x0)/gamma;
switch diff
  case -1
    yabs = 1/2*(1+erf(sqrt(2)*k));
    ydisp = NaN;
  case 0
    yabs = pre/gamma*exp(-2*k.^2);
    if (nargout>1)
      %dawson = @(x) real(sqrt(pi)/2i*(faddeeva(x)-exp(-x.^2))); % general x
      dawson = @(x) sqrt(pi)/2*imag(faddeeva(x)); % for real x
      ydisp = pre/gamma*2/sqrt(pi)*dawson(sqrt(2)*k);
    end
  case 1
    yabs = -pre*4/gamma^2*k.*exp(-2*k.^2);
    ydisp = NaN;
  case 2
    yabs = pre*4/gamma^3*(4*k.^2-1).*exp(-2*k.^2);
    ydisp = NaN;
end

return
