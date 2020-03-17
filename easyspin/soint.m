% soint  Spin-orbit interaction Hamiltonian 
%
%   H = soint(SpinSystem)
%   H = soint(SpinSystem,Spins)
%   H = soint(SpinSystem,Spins,'sparse')
%
%   Returns the spin orbit interaction (SOI)
%   Hamiltonian, in MHz.
%
%   Input:
%   - SpinSystem: Spin system structure. SOI
%       parameters are the spin-orbit coupling in soc 
%       and the orbital reduction facotr in orf.
%   - Spins: If given, specifies spins
%       for which the SOI should be computed. If
%       absent, all electrons are included.
%   - 'sparse': If given, the matrix is returned in sparse format.
%
%   Output:
%   - H: Hamiltonian matrix containing the SOI for
%        spins specified in eSpins.

function H = soint(System,Spins,opt)

if nargin==0, help(mfilename); return; end

if nargin<1 || nargin>3, error('Wrong number of input arguments!'); end
if nargout>1, error('Too many output arguments.'); end

if nargin<3, opt = ''; end
if nargin<2, Spins = []; end
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

if isempty(Spins), Spins = 1:System.nElectrons; end

% Some error checking on the second input argument
if numel(Spins)<1
  error('Spins (2nd argument) must contain at least 1 values!'); 
end
if any(Spins<1) || any(Spins>System.nElectrons)
  error('Spins (2nd argument) must contain values between 1 and %d!',System.nElectrons);
end
if numel(unique(Spins))~=numel(Spins)
  error('Spins (2nd argument) contains double entries!');
end


% Spin-orbit coupling terms soc*S(k)*L(k)
%----------------------------------------------------------------
shift = System.nElectrons + System.nNuclei;
for k = 1:length(Spins)
  
  iSpin = Spins(k);
  soc_ = System.soc(iSpin,:);
  if ~any(soc_), continue; end
  
  % Build S*L operator matrix
  SL = sparse(System.nStates,System.nStates);
  for c1 = 1:3
    SL = SL + sop(System,[iSpin c1; iSpin+shift c1],'sparse');
  end
  
  % Add SOC terms to Hamiltonian
  for order = find(soc_)
    H = H + soc_(order) * (System.orf(iSpin)*SL)^order;
  end
  
end

if ~sparseResult
  H = full(H); % sparse -> full
end
