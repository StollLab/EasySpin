% hsdim  Dimension on spin state space
%
%   nStates = hsdim(SpinSystem)
%
%   Returns the Hilbert space dimension of spin system 'SpinSystem', which can
%   be either a spin system structure or a vector containing spin quantum
%   numbers.
%
%   All electron spins, nuclear spins (including classes of equivalent ones),
%   and orbital angular momenta are included in the calculation.

function nStates = hsdim(SpinSystem)

if nargin==0, help(mfilename); return; end

if isstruct(SpinSystem)
  [Sys,err] = validatespinsys(SpinSystem);
  error(err);
  S = Sys.S;
  I = Sys.I;
  L = Sys.L;
  if isfield(SpinSystem,'n')
    nEquivI = SpinSystem.n;
  else
    nEquivI = [];
  end
else
  S = SpinSystem;
  L = [];
  I = [];
  nEquivI = [];
end

if isempty(nEquivI)
  J = [S I L];
  nStates = prod(2*J+1);
else
  nStates = prod(2*S+1) * prod((2*I+1).^nEquivI) * prod(2*L+1);
end

return
