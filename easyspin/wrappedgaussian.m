% wrappedgaussian  "Wrapped" Gaussian function
%
%   y = wrappedgaussian(x,domain,x0,fwhm)
%
%   Computes area-normalized Gaussian function that is "wrapped" to satisfy
%   periodic boundary conditions
%
%   Input:
%   - x:        abscissa vector
%   - x0:       center of the function
%   - fwhm:     full width at half height
%   - domain:   minimum and maximum x-values of domain to be wrapped, e.g.
%               use [-pi, pi] to define the function's domain to be 
%               -pi <= x <= pi
%
%   Output:
%   - y:   function values for abscissa x

% Reference:
%  See https://en.wikipedia.org/wiki/Wrapped_normal_distribution

function y = wrappedgaussian(x,x0,fwhm,domain)

if (nargin==0), help(mfilename); return; end

if (nargin<4), domain = [-pi,pi]; end

if any(fwhm<=0) || any(~isreal(fwhm))
  error('fwhm must be positive and real!');
end

if numel(fwhm)>1
  error('fwhm must contain 1 element!');
end

if ~isrow(domain)||numel(domain)~=2||any(~isreal(domain))
  error('domain must be a 2-element row vector of real numbers.')
end

if domain(1)>=domain(2)
  error('domain elements must be given in increasing order.')
end

N = 5;
shift = domain(2)-domain(1);

k = (-5:5).';  % column of k's to be used for summation
xgrid = x - x0 + shift*k;  % grid of argument domains to be evaluated in exponential
y = sum(gaussian(xgrid, 0, fwhm),1);

end



