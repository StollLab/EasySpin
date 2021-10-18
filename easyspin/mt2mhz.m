% mt2mhz   Conversion from mT to MHz
%
%   v_MHz = mt2mhz(B_mT)
%   v_MHz = mt2mhz(B_mT,g)
%
%   Converts the magnetic field value in B_mT, assumed to be
%   in units of mT (Millitesla), to MHz (Megahertz), returning
%   the result in v_MHz.
%
%   For the conversion, the g factor given as second
%   parameter is used. If it is not given, the g factor
%   of the free electron (gfree) is used.
%
%   B_mT can be a vector of values. In this case, g
%   can be a scalar or a vector of the same size as B_mT.

function x_MHz = mt2mhz(x_mT,g)

if (nargin<1), x_mT = 1; end
if (nargin<2), g = gfree; end

x_MHz = x_mT.*g*(1e-3*bmagn/planck/1e6);

return
