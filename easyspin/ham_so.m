% ham_so  Spin-orbit interaction Hamiltonian 
%
%   H = ham_so(SpinSystem)
%   H = ham_so(SpinSystem,eSpins)
%   H = ham_so(SpinSystem,eSpins,'sparse')
%
%   Returns the spin orbit interaction (SOI) Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. The spin-orbit coupling is
%       in SpinSystem.soc.
%   - Spins: If given, specifies electron spins for which the SOI should be
%       computed. If absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.
%
%   Output:
%   - H: Hamiltonian matrix containing the SOI for spins specified in eSpins.

function H = ham_so(System,eSpins,opt)

if nargin==0, help(mfilename); return; end

if nargin<1 || nargin>3, error('Wrong number of input arguments!'); end
if nargout>1, error('Too many output arguments.'); end

if nargin<3, opt = ''; end
if nargin<2, eSpins = []; end
if ~ischar(opt)
  error('Third input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

[System,err] = validatespinsys(System);
error(err);

% Special cases: no orbital angular momemta L, or .soc not given or all zero
H = sparse(System.nStates,System.nStates);
if System.nL==0 || all(System.soc(:)==0)
  if ~sparseResult
    H = full(H); % sparse -> full
  end
  return
end

if isempty(eSpins), eSpins = 1:System.nElectrons; end

% Some error checking on the second input argument
if numel(eSpins)<1
  error('eSpins (2nd argument) must contain at least 1 values!'); 
end
if any(eSpins<1) || any(eSpins>System.nElectrons)
  error('eSpins (2nd argument) must contain values between 1 and %d!',System.nElectrons);
end
if numel(unique(eSpins))~=numel(eSpins)
  error('eSpins (2nd argument) contains double entries!');
end


% Spin-orbit coupling terms soc*S(k)*L(k)
%--------------------------------------------------------------------------
shift = System.nElectrons + System.nNuclei;
for k = 1:length(eSpins)
  
  iS = eSpins(k);
  iL = iS+shift;
  soc_ = System.soc(iS,:);
  if ~any(soc_), continue; end

  % Skip if orbital angular momentum is zero
  if System.L(iS)==0, continue; end
  
  % Build S*L operator matrix
  SL = sparse(System.nStates,System.nStates);
  for c1 = 1:3
    SL = SL + sop(System,[iS c1; iL c1],'sparse');
  end
  
  % Add SOC terms to Hamiltonian
  for order = find(soc_)
    H = H + soc_(order) * SL^order;
  end
  
end

if ~sparseResult
  H = full(H); % sparse -> full
end

end
