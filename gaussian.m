% gaussian  Gaussian line shape 
%
%   y = gaussian(x,x0,fwhm)
%   y = gaussian(x,x0,fwhm,diff)
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
%   - y:    Vector of function values for arguments x

function y = gaussian(x,x0,fwhm,diff)

if (nargin==0)
  if (nargout==0)
    help(mfilename);
    return;
  else
    x = linspace(-1,1,1001);
    y = gaussian(x,0,0.2);
    return;
  end
end

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

% prefactor is 1/sqrt(2*log(2))
gamma = 0.849321800288*fwhm; % distance from x0 to inflexion point
pre = 0.797884560803;  % sqrt(2/pi)
k = (x-x0)/gamma;
switch diff
  case -1
    y = 1/2*(1+erf(sqrt(2)*k));
  case 0
    y = pre/gamma*exp(-2*k.^2);
  case 1
    y = -pre*4/gamma^2*k.*exp(-2*k.^2);
  case 2
    y = pre*4/gamma^3*(4*k.^2-1).*exp(-2*k.^2);
end

return