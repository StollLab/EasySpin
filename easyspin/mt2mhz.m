% mt2mhz   Conversion from mT to MHz
%
%   value_MHz = mt2mhz(value_mT)
%   value_MHz = mt2mhz(value_mT,g)
%
%   Converts the magnetic field value in value_mT, assumed to be
%   in units of mT (millitesla), to MHz (megahertz), returning
%   the result in value_MHz.
%
%   For the conversion, the g factor given as second
%   parameter is used. If it is not given, the g factor
%   of the free electron (gfree) is used.
%
%   value_mT can be a vector of values. In this case, g
%   can be a scalar or a vector of the same size as value_mT.

function value_MHz = mt2mhz(value_mT,g)

if nargin==0, help(mfilename); return; end

if nargin<2, g = gfree; end

value_MHz = unitconvert(value_mT,"mT->MHz",g);

end
