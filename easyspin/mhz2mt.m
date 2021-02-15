% mhz2mt   Conversion from MHz to mT
%
%   x_mT = mhz2mt(x_MHz)
%   x_mT = mhz2mt(x_MHz,g)
%
%   Converts the value in x_MHz, assumed to be in units
%   of MHz (megahertz) to mT (millitesla), returning
%   the result in x_mT.
%
%   For the conversion, the g factor given as second
%   parameter is used. If it is not given, the g factor
%   of the free electron (gfree) is used.
%
%   x_MHz can be a vector of values. In this case, g
%   can be a scalar or a vector of the same size as x_MHz.

function x_mT = mhz2mt(x_MHz,g)

if (nargin<1), x_MHz = 1; end
if (nargin<2), g = gfree; end

if ~isnumeric(g)
  error('Second input (g) must be numeric.');
end

x_mT = x_MHz./g*(1e6*planck/bmagn/1e-3);

return
