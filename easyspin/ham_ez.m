% ham_ez  Electron Zeeman interaction Hamiltonian 
%
%   H = ham_ez(SpinSystem, B)
%   H = ham_ez(SpinSystem, B, eSpins)
%   H = ham_ez(SpinSystem, B, eSpins, 'sparse')
%   [Gx, Gy, Gz] = ham_ez(SpinSystem)
%   [Gx, Gy, Gz] = ham_ez(SpinSystem, eSpins)
%   [Gx, Gy, Gz] = ham_ez(SpinSystem, eSpins, 'sparse')
%
%   Returns the electron Zeeman interaction Hamiltonian matrix for
%   the electron spins 'eSpins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - eSpins: Vector of indices for electron spins to include. If eSpins is
%     omitted, all electron spins are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - Gx, Gy, Gz: Components of the Zeeman interaction Hamiltonian matrix
%     for the selected electron spins as defined by Gi=d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field in the molecular frame. Units are MHz/mT = GHz/T.
%     To get the full electron Zeeman Hamiltonian, use H = Gx*B(1)+Gy*B(2)+Gz*B(3),
%     where B is the magnetic field vector in mT.
%   - H: Electron Zeeman Hamiltonian matrix.

function varargout = ham_ez(SpinSystem,varargin)

if nargin==0, help(mfilename); return; end

if nargout~=1 && nargout~=3, error('Wrong number of output arguments!'); end
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
GxM = sparse(nStates,nStates);
GyM = sparse(nStates,nStates);
GzM = sparse(nStates,nStates);

% Calculate prefactors
pre = +bmagn/planck*Sys.g; % Hz/T
pre = pre/1e9;  % Hz/T -> GHz/T = MHz/mT

% Loop over all electron spins
for i = eSpins

  if Sys.fullg
    g = pre((i-1)*3+(1:3),:);
  else
    g = diag(pre(i,:));
  end
  % Transform g matrix to molecular frame
  R_M2g = erot(Sys.gFrame(i,:));  % mol frame -> g frame
  R_g2M = R_M2g.';  % g frame -> mol frame
  g = R_g2M*g*R_g2M.';
  % Build electon Zeeman Hamiltonian in MHz/mT
  for k = 1:3
    Sk = sop(spins,[i,k],'sparse');
    GxM = GxM + g(1,k)*Sk;
    GyM = GyM + g(2,k)*Sk;
    GzM = GzM + g(3,k)*Sk;
  end

end

if isempty(B0)
  if ~useSparseMatrices
    GxM = full(GxM);
    GyM = full(GyM);
    GzM = full(GzM);
  end
  varargout = {GxM, GyM, GzM};
else
  H = GxM*B0(1) + GyM*B0(2) + GzM*B0(3);
  if ~useSparseMatrices
    H = full(H);
  end
  varargout = {H};
end

end
