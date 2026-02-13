% ham_ez  Electron Zeeman interaction Hamiltonian 
%
%   H = ham_ez(SpinSystem, B)
%   H = ham_ez(SpinSystem, B, eSpins)
%   H = ham_ez(SpinSystem, B, eSpins, 'sparse')
%   [mux, muy, muz] = ham_ez(SpinSystem)
%   [mux, muy, muz] = ham_ez(SpinSystem, eSpins)
%   [mux, muy, muz] = ham_ez(SpinSystem, eSpins, 'sparse')
%   [mux, muy, muz, dmudgx, dmudgy, dmudgz] = ham_ez(SpinSystem)
%   [mux, muy, muz, dmudgx, dmudgy, dmudgz] = ham_ez(SpinSystem, eSpins)
%   [mux, muy, muz, dmudgx, dmudgy, dmudgz] = ham_ez(SpinSystem, eSpins, 'sparse')
%
%   Returns the electron Zeeman interaction Hamiltonian matrix for
%   the electron spins 'eSpins' of the spin system 'SpinSystem'. It can
%   return also the magnetic moment and its derivatives with respect to gx,
%   gy, gz in the molecular frame.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - eSpins: Vector of indices for electron spins to include. If eSpins is
%     omitted, all electron spins are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - mux, muy, muz: Components of the magnetic dipole moment operator
%     for the selected electron spins as defined by mui=-d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full electron Zeeman Hamiltonian, use
%     H = -(mux*B(1)+muy*B(2)+muz*B(3)), where B is the magnetic field vector
%     in mT.
%   - dmudgx, dmudgy, dmudgz, derivatives of the magnetic dipole moment
%     with respect to dmudgi= dmu/dgi where i=x,y,z stored as a cell index as {spin,i}
%   - H: Electron Zeeman Hamiltonian matrix.

function varargout = ham_ez(SpinSystem,varargin)

if nargin==0, help(mfilename); return; end

if nargout~=1 && nargout~=3 && nargout~=6
  error('Wrong number of output arguments!');
end
singleOutput = nargout==1;

if singleOutput
  if nargin<1 || nargin>4
    error('Wrong number of input arguments!');
  end
  if nargin<2
    error('Field vector (second input, in mT) is missing.')
  else
    B0 = varargin{1};
  end
  if nargin<3, eSpins = []; else, eSpins = varargin{2}; end
  if nargin<4, opt = ''; else, opt = varargin{3}; end
else
  B0 = [];
  if nargin<2, eSpins = []; else, eSpins = varargin{1}; end
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

% No 'eSpins' specified -> use all
if isempty(eSpins), eSpins = 1:nElectrons; end

% Validate magnetic field if given
if ~isempty(B0)
  if numel(B0)~=3
    error('Magnetic field vector (2nd input) must be a 3-element array.');
  end
end

% Validate third argument (eSpins)
if any(eSpins<1) || any(eSpins>nElectrons)
  error('Electron spin indices (2nd input argument) invalid!');
end

% Initialize Zeeman interaction component matrices to zero
nStates = Sys.nStates;
muxM = sparse(nStates,nStates);
muyM = sparse(nStates,nStates);
muzM = sparse(nStates,nStates);

% Calculate prefactors
pre = -bmagn/planck*Sys.g; % Hz/T
pre = pre/1e9;  % Hz/T -> GHz/T = MHz/mT
preDer = -bmagn/planck/1e9; 

% Loop over all electron spins
for eSp = eSpins

  % Get full g matrix
  if Sys.fullg
    g = pre((eSp-1)*3+(1:3),:);
  else
    g = diag(pre(eSp,:));
  end

  % Transform g matrix to molecular frame
  ang = Sys.gFrame(eSp,:);
  if any(ang)
    R_M2g = erot(ang);  % mol frame -> g frame
    R_g2M = R_M2g.';  % g frame -> mol frame
    g = R_g2M*g*R_g2M.';
  else
    R_g2M=eye(3);
  end

  % preparing the derivatives of the magnetic moment
  for index=1:3
    dmuMgx{index} = sparse(nStates,nStates);
    dmuMgy{index} = sparse(nStates,nStates);
    dmuMgz{index} = sparse(nStates,nStates);
  end

  % preparing the derivatives coefficient
  dgxM = preDer*R_g2M(:,1)*R_g2M(:,1).';  % rotate derivative wrt gx to molecular frame (and scale with preDer)
  dgyM = preDer*R_g2M(:,2)*R_g2M(:,2).';  % rotate derivative wrt gy to molecular frame (and scale with preDer)
  dgzM = preDer*R_g2M(:,3)*R_g2M(:,3).';  % rotate derivative wrt gz to molecular frame (and scale with preDer)

  % Build magnetic dipole moment components in MHz/mT
  for k = 1:3
    if ~useSparseMatrices
      Sk = sop(spins,[eSp,k]);
    else
      Sk = sop(spins,[eSp,k],'sparse');
    end
    muxM = muxM + g(1,k)*Sk;
    muyM = muyM + g(2,k)*Sk;
    muzM = muzM + g(3,k)*Sk;

    for index=1:3
      dmuMgx{index} = dmuMgx{index} + dgxM(index,k)*Sk;
      dmuMgy{index} = dmuMgy{index} + dgyM(index,k)*Sk;
      dmuMgz{index} = dmuMgz{index} + dgzM(index,k)*Sk;
    end
  end

    dmuMdgx{eSp} = {dmuMgx{1}, dmuMgx{2}, dmuMgx{3}};
    dmuMdgy{eSp} = {dmuMgy{1}, dmuMgy{2}, dmuMgy{3}};
    dmuMdgz{eSp} = {dmuMgz{1}, dmuMgz{2}, dmuMgz{3}};
end

if isempty(B0)
  % Return magnetic dipole moment components and derivatives
  varargout = {muxM, muyM, muzM, dmuMdgx, dmuMdgy, dmuMdgz};
else
  % Return Zeeman Hamiltonian
  H = -(muxM*B0(1) + muyM*B0(2) + muzM*B0(3));
  varargout = {H};
end

end
