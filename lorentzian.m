% lorentzian  Lorentzian line shape 
%
%   ya = lorentzian(x,x0,fwhm)
%   ya = lorentzian(x,x0,fwhm,diff)
%   [ya,yd] = ...
%
%   Returns area-normalized Lorentzian absorption and
%   dispersion line shapes or their derivatives.
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
%   - ya:    absorption function values for abscissa x
%   - yd:    dispersion function values for abscissa x

function [yabs,ydisp] = lorentzian(x,x0,fwhm,diff)

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

% Compute Lorentzian lineshape
%------------------------------------------------------------------
% gamma = distance from center to inflexion point
gamma = fwhm/sqrt(3);
pre = 2/pi/sqrt(3);
k = (x-x0)/gamma;
switch diff
  case -1
    yabs = atan(2/sqrt(3)*k)/pi + 1/2;
    ydisp = NaN;
  case 0
    yabs = pre/gamma./(1+4/3*k.^2);
    ydisp = -pre^2*pi/gamma*k./(1+4/3*k.^2);
  case 1
    yabs = -8/3*pre/gamma^2*k./(1+4/3*k.^2).^2;
    ydisp = -pre^2*pi/gamma^2*(1-4/3*k.^2)./(1+4/3*k.^2).^2;
  case 2
    yabs = 8/3*pre/gamma^3*(4*k.^2-1)./(1+4/3*k.^2).^3;
    ydisp = -pre^2*pi/gamma^3*2*4/3*k.*(4/3*k.^2-3)./(1+4/3*k.^2).^3;
end

return
