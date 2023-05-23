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
%            2 second and so on, -1 the integral with -infinity
%            as lower limit. 0 is the default.
%   - phase: phase rotation, mixes absorption and dispersion.
%            phase=pi/2 puts dispersion signal into ya
%
%   Output:
%   - ya:    absorption function values for abscissa x
%   - yd:    dispersion function values for abscissa x

function [yabs,ydisp] = lorentzian(x,x0,fwhm,diff,phase)

if nargin==0, help(mfilename); return; end

if nargin<4, diff = 0; end
if nargin<5, phase = 0; end

if any(fwhm<=0) || any(~isreal(fwhm))
  error('fwhm must be positive and real!');
end

if numel(fwhm)>1
  error('fwhm must contain 1 element!');
end

if any(diff<-1) || any(diff>2)
  error('Cannot compute Lorentzian lineshape derivative %d.',diff);
end

if numel(phase)>1
  error('phase must contain 1 element.');
elseif ~isreal(phase)
  error('phase must be real.');
end


% Compute Lorentzian lineshape
%-------------------------------------------------------------------------------
gamma = fwhm/sqrt(3);  % distance from center to inflection point
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
  yabs_ =  yabs*cos(phase) + ydisp*sin(phase);
  ydisp = -yabs*sin(phase) + ydisp*cos(phase);
  yabs = yabs_;
end

end
