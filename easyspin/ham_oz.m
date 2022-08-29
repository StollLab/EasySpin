% ham_ez  Orbital Zeeman interaction Hamiltonian 
%
%   H = ham_oz(SpinSystem, B)
%   H = ham_oz(SpinSystem, B, oam)
%   H = ham_oz(SpinSystem, B, oam, 'sparse')
%   [Zx, Zy, Zz] = ham_oz(SpinSystem)
%   [Zx, Zy, Zz] = ham_oz(SpinSystem, oam)
%   [Zx, Zy, Zz] = ham_oz(SpinSystem, oam, 'sparse')
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
%   - Zx, Zy, Zz: Components of the Zeeman interaction Hamiltonian matrix
%     for the selected orbital angular momenta as defined by Hi=d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full Hamiltonian, use H = Zx*B(1)+Zy*B(2)+Zz*B(3), where
%     B is the magnetic field in mT.
%   - H: Orbital Zeeman Hamiltonian matrix.

function varargout = ham_oz(SpinSystem,varargin)

if nargin==0, help(mfilename); return; end

if nargout==2 || nargout>3, error('Wrong number of output arguments!'); end
singleOutput = nargout<2;

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
ZxM = sparse(nStates,nStates);
ZyM = sparse(nStates,nStates);
ZzM = sparse(nStates,nStates);

% Calculate prefactor
pre = +bmagn/planck*Sys.gL; % Hz/T
pre = pre/1e9;  % Hz/T -> GHz/T = MHz/mT

% Loop over all orbital angular momenta
for i = oam

  % Orbital angular momenta, isotropic
  % Build orbital Zeeman Hamiltonian in MHz/mT
  iSpin = nElectrons + nNuclei + i;
  ZxM = ZxM + pre(i)*sop(spins,[iSpin,1],'sparse');
  ZyM = ZyM + pre(i)*sop(spins,[iSpin,2],'sparse');
  ZzM = ZzM + pre(i)*sop(spins,[iSpin,3],'sparse');

end

if isempty(B0)
  if ~useSparseMatrices
    ZxM = full(ZxM);
    ZyM = full(ZyM);
    ZzM = full(ZzM);
  end
  varargout = {ZxM, ZyM, ZzM};
else
  H = ZxM*B0(1) + ZyM*B0(2) + ZzM*B0(3);
  if ~useSparseMatrices
    H = full(H);
  end
  varargout = {H};
end

end
