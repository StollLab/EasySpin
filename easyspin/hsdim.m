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
  [Sys,err] = validatespinsys(SpinSystem);
  error(err);
  eSpins = Sys.S;
  nSpins = Sys.I;
  if isfield(SpinSystem,'n')
    nEquiv = SpinSystem.n;
  else
    nEquiv = [];
  end
else
  eSpins = SpinSystem;
  nSpins = [];
  nEquiv = [];
end

if isempty(nEquiv)
  nStates = prod(2*[eSpins nSpins]+1);
else
  nStates = prod(2*eSpins+1)*prod((2*nSpins+1).*nEquiv);
end

return
