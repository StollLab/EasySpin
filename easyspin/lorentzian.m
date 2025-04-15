% lorentzian  Lorentzian line shape 
%
%   ya = lorentzian(x,x0,fwhm)
%   ya = lorentzian(x,x0,fwhm,diff)
%   ya = lorentzian(x,x0,fwhm,diff,phase)
%   [ya,yd] = lorentzian(...)
%
%   Computes area-normalized Lorentzian absorption and dispersion
%   line shapes or their derivatives.
%
%   Input:
%   - x:     abscissa vector
%   - x0:    center of the lineshape function
%   - fwhm:  full width at half height
%   - diff:  derivative. 0 is no derivative, 1 first,
%            2 second derivative, -1 the integral with -infinity
%            as lower limit. 0 is the default.
%   - phase: phase rotation, mixes absorption and dispersion.
%            phase=pi/2 puts dispersion signal into ya
%
%   Output:
%   - ya:    absorption function values for abscissa x
%   - yd:    dispersion function values for abscissa x

function [y1,y2] = lorentzian(x,x0,fwhm,diff,phase)

if nargin==0, help(mfilename); return; end

% Check number of input arguments
%-------------------------------------------------------------------------------
if nargin<3
  error('lorentzian() needs at least three inputs: x, x0, and fwhm.');
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

if ~isscalar(diff) || diff<-1 || diff>2 || mod(diff,1)
  error('4th input (diff, lineshape derivative) must be -1, 0, 1, or 2.');
end

if ~isscalar(phase) || ~isreal(phase)
  error('5th input (phase) must be a real number.');
end

% Compute Lorentzian lineshape
%-------------------------------------------------------------------------------
% gamma = distance from center to inflection point
gamma = fwhm/sqrt(3);
pre = 2/pi/sqrt(3);
z = (x-x0)/gamma;
switch diff
  case -1
    yabs = atan(2/sqrt(3)*z)/pi + 1/2;
    ydisp = (1/2/pi)*log(3+4*z.^2);
  case 0
    yabs = pre/gamma./(1+4/3*z.^2);
    ydisp = pre^2*pi/gamma*z./(1+4/3*z.^2);
  case 1
    yabs = -8/3*pre/gamma^2*z./(1+4/3*z.^2).^2;
    ydisp = pre^2*pi/gamma^2*(1-4/3*z.^2)./(1+4/3*z.^2).^2;
  case 2
    yabs = 8/3*pre/gamma^3*(4*z.^2-1)./(1+4/3*z.^2).^3;
    ydisp = pre^2*pi/gamma^3*2*4/3*z.*(4/3*z.^2-3)./(1+4/3*z.^2).^3;
end

% Phase rotation
%-------------------------------------------------------------------------------
if phase~=0
  y1 =  yabs*cos(phase) + ydisp*sin(phase);
  y2 = -yabs*sin(phase) + ydisp*cos(phase);
else
  y1 = yabs;
  y2 = ydisp;
end

end
