% ham_nz  Nuclear Zeeman interaction Hamiltonian 
%
%   H = ham_nz(SpinSystem, B)
%   H = ham_nz(SpinSystem, B, nSpins)
%   H = ham_nz(SpinSystem, B, nSpins, 'sparse')
%   [mux, muy, muz] = ham_nz(SpinSystem)
%   [mux, muy, muz] = ham_nz(SpinSystem, nSpins)
%   [mux, muy, muz] = ham_nz(SpinSystem, nSpins, 'sparse')
%
%   Returns the nuclear Zeeman interaction Hamiltonian matrix for
%   the nuclear spins 'nSpins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - nSpins: Vector of indices for nuclear spins to include. If nSpins is
%     omitted, all nuclear spins are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - mux, muy, muz: Components of the magnetic dipole moment operator
%     for the selected nuclear spins as defined by mui=-d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full nuclear Zeeman Hamiltonian, use
%     H = -(mux*B(1)+muy*B(2)+muz*B(3)), where B is the magnetic field vector in mT.
%   - H: Nuclear Zeeman Hamiltonian matrix.

function varargout = ham_nz(SpinSystem,varargin)

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
  if nargin<3, nSpins = []; else, nSpins = varargin{2}; end
  if nargin<4, opt = ''; else, opt = varargin{3}; end
else
  B0 = [];
  if nargin<2, nSpins = []; else, nSpins = varargin{1}; end
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
nNuclei = Sys.nNuclei;

% No 'eSpins' specified -> use all
if isempty(nSpins), nSpins = 1:nNuclei; end

% Validate magnetic field if given
if ~isempty(B0)
  if numel(B0)~=3
    error('Magnetic field vector (2nd input) must be a 3-element array.');
  end
end

% Validate third argument (nSpins)
if any(nSpins<1) || any(nSpins>nNuclei)
  error('Nuclear spin indices (2nd input argument) invalid!');
end

% Initialize Zeeman interaction component matrices to zero
nStates = Sys.nStates;
muxM = sparse(nStates,nStates);
muyM = sparse(nStates,nStates);
muzM = sparse(nStates,nStates);

% Calculate prefactors
pre = +nmagn/planck*Sys.gn.*Sys.gnscale; % Hz/T
pre = pre/1e9;  % Hz/T -> GHz/T = MHz/mT

% Loop over all nuclear spins
nElectrons = Sys.nElectrons;
for iNuc = nSpins

  % Nuclei, with isotropic gn and chemical shielding (CS) tensor sigma
  if Sys.fullsigma
    sigma = Sys.sigma((iNuc-1)*3+(1:3),:);
  else
    sigma = diag(Sys.sigma(iNuc,:));
  end

  % Transform CS tensor to molecular frame
  ang = Sys.sigmaFrame(iNuc,:);
  if any(ang)
    R_M2CS = erot(ang);  % mol frame -> CS frame
    R_CS2M = R_M2CS.';   % CS frame -> mol frame
    sigma = R_CS2M*sigma*R_CS2M.';
  end

  % Build nuclear Zeeman Hamiltonian in MHz/mT
  for k = 1:3
    Ik = sop(spins,[nElectrons+iNuc,k],'sparse');
    muxM = muxM + pre(iNuc)*sigma(1,k)*Ik;
    muyM = muyM + pre(iNuc)*sigma(2,k)*Ik;
    muzM = muzM + pre(iNuc)*sigma(3,k)*Ik;
  end

end

if isempty(B0)
  if ~useSparseMatrices
    muxM = full(muxM);
    muyM = full(muyM);
    muzM = full(muzM);
  end
  varargout = {muxM, muyM, muzM};
else
  H = muxM*B0(1) + muyM*B0(2) + muzM*B0(3);
  if ~useSparseMatrices
    H = full(H);
  end
  varargout = {H};
end

end
