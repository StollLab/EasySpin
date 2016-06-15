% hsdim  Dimension on spin state space 
%
%   nStates = hsdim(SpinSystem)
%
%   Returns the Hilbert space
%   dimension of spin system 'SpinSystem',
%   which can be either a spin system
%   structure or a vector containing
%   the spin quantum numbers

function nStates = hsdim(SpinSystem)

if (nargin==0), help(mfilename); return; end

if isstruct(SpinSystem)
  SpinSystem = spinvec(SpinSystem);
end

nStates = prod(2*SpinSystem + 1);

return
