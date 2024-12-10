% ham_ez  Orbital Zeeman interaction Hamiltonian 
%
%   H = ham_oz(SpinSystem, B)
%   H = ham_oz(SpinSystem, B, oam)
%   H = ham_oz(SpinSystem, B, oam, 'sparse')
%   [mux, muy, muz] = ham_oz(SpinSystem)
%   [mux, muy, muz] = ham_oz(SpinSystem, oam)
%   [mux, muy, muz] = ham_oz(SpinSystem, oam, 'sparse')
%
%   Returns the orbital Zeeman interaction Hamiltonian matrix for
%   the orbital angular momenta 'oam' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - oam: Vector of indices for orbital angular momenta to include. If oam
%     is omitted, all orbital angular momenta are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - mux, muy, muz: Components of the magnetic dipole moment operator
%     for the selected orbital angular momenta as defined by mui=d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full orbital Zeeman Hamiltonian, use
%     H = -(mux*B(1)+muy*B(2)+muz*B(3)), where B is the magnetic field vector in mT.
%   - H: Orbital Zeeman Hamiltonian matrix.

function varargout = ham_oz(SpinSystem,varargin)

if nargin==0, help(mfilename); return; end

if nargout~=1 && nargout~=3, error('Wrong number of output arguments!'); end
singleOutput = nargout==1;

if singleOutput
  if nargin<1 || nargin>4, error('Wrong number of input arguments!'); end
  if nargin<2
    error('Field vector (second input, in mT) is missing.')
  else
    B0 = varargin{1};
  end
  if nargin<3, oam = []; else, oam = varargin{2}; end
  if nargin<4, opt = ''; else, opt = varargin{3}; end
else
  B0 = [];
  if nargin<2, oam = []; else, oam = varargin{1}; end
  if nargin<3, opt = ''; else, opt = varargin{2}; end
end

if ~ischar(opt)
  error('Last input must be a string, ''sparse''.');
end
useSparseMatrices = strcmp(opt,'sparse');

% Validate spin system
[Sys,err] = validatespinsys(SpinSystem);
error(err);

% Vector of all spin quantum numbers
spins = Sys.Spins;
nElectrons = Sys.nElectrons;
nNuclei = Sys.nNuclei;
nOam = Sys.nL;

% No 'oam' specified -> use all
if isempty(oam), oam = 1:nOam; end

% Validate magnetic field if given
if ~isempty(B0)
  if numel(B0)~=3
    error('Magnetic field vector (2nd input) must be a 3-element array.');
  end
end

% Validate third argument (oam)
if any(oam<1) || any(oam>nElectrons)
  error('Orbital angular momentum indices (2nd input argument) invalid!');
end

% Initialize Zeeman interaction component matrices to zero
nStates = Sys.nStates;
muxM = sparse(nStates,nStates);
muyM = sparse(nStates,nStates);
muzM = sparse(nStates,nStates);

% Calculate prefactor
pre = -bmagn/planck*Sys.gL; % Hz/T
pre = pre/1e9;  % Hz/T -> GHz/T = MHz/mT

% Loop over all orbital angular momenta
for i = oam

  % Skip if orbital angular momenutm is zero
  if Sys.L(i)==0; continue; end

  % Orbital angular momenta, isotropic
  % Build orbital Zeeman Hamiltonian in MHz/mT
  iSpin = nElectrons + nNuclei + i;
  muxM = muxM + pre(i)*sop(spins,[iSpin,1],'sparse');
  muyM = muyM + pre(i)*sop(spins,[iSpin,2],'sparse');
  muzM = muzM + pre(i)*sop(spins,[iSpin,3],'sparse');

end

if isempty(B0)
  if ~useSparseMatrices
    muxM = full(muxM);
    muyM = full(muyM);
    muzM = full(muzM);
  end
  varargout = {muxM, muyM, muzM};
else
  H = -(muxM*B0(1) + muyM*B0(2) + muzM*B0(3));
  if ~useSparseMatrices
    H = full(H);
  end
  varargout = {H};
end

end
