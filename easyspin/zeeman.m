% zeeman  Zeeman interaction Hamiltonian 
%
%   H = zeeman(SpinSystem, B)
%   H = zeeman(SpinSystem, B, idx)
%   H = zeeman(SpinSystem, B, idx, 'sparse')
%   [Zx, Zy, Zz] = zeeman(SpinSystem)
%   [Zx, Zy, Zz] = zeeman(SpinSystem, idx)
%   [Zx, Zy, Zz] = zeeman(SpinSystem, idx, 'sparse')
%
%   Returns the Zeeman interaction Hamiltonian for
%   the spins 'Spins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - idx: Vector of indices for spins/orbital angular momenta to include.
%     For one electron spin: 1 is the electron, >=2 are the nuclei. For two
%     electron spins: 1 and 2 electrons, >=3 nuclei, etc. If also orbital
%     angular momenta are defined, they follow after the nuclei. If idx is
%     omitted, all spins and orbital angular momenta are included.
%     Two electron spins, three nuclei and 2 orbital angular momenta: 1 and 2
%     are electrons, 3,4, and 5 are nuclei, 6 and 7 are orbital angular momenta. 
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - Zx, Zy, Zz: components of the Zeeman interaction Hamiltonian
%     for the selected spins and orbital angular momenta as defined by Hi=d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field. Units are MHz/mT = GHz/T. To get the
%     full Hamiltonian, use H = Zx*B(1)+Zy*B(2)+Zz*B(3), where B is
%     the magnetic field in mT.
%   - H: the Hamiltonian of the Zeeman interaction.

function varargout = zeeman(SpinSystem,varargin)

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
  if nargin<3, idx = []; else, idx = varargin{2}; end
  if nargin<4, opt = ''; else, opt = varargin{3}; end
else
  B0 = [];
  if nargin<2, idx = []; else, idx = varargin{1}; end
  if nargin<3, opt = ''; else, opt = varargin{2}; end
end

if ~ischar(opt)
  error('Last input must be a string, ''sparse''.');
end
useSparseMatrices = strcmp(opt,'sparse');

% Validate spin system
[Sys,err] = validatespinsys(SpinSystem);
error(err);

% Vector of angular momemtum quantum numbers (e-spin S, n-spin I, e-orbital L)
QuantumNumbers = Sys.Spins;

% No 'idx' specified -> use all
if isempty(idx), idx = 1:numel(QuantumNumbers); end

% Validate magnetic field if given
if ~isempty(B0)
  if numel(B0)~=3
    error('Magnetic field vector (2nd input) must be a 3-element array.');
  end
end

% Validate third argument (idx)
if any(idx<1) || any(idx>length(QuantumNumbers))
  error('Indices (2nd input argument) invalid!');
end

% Get number of electron spins, nuclear spins, and orrbital angular momenta
nElectrons = Sys.nElectrons;
nNuclei = Sys.nNuclei;
nStates = Sys.nStates;

% Initialize Zeeman interaction component matrices to zero
ZxM = sparse(nStates,nStates);
ZyM = sparse(nStates,nStates);
ZzM = sparse(nStates,nStates);

% Complete prefactors, in MHz/mT
elFactor = +bmagn/(planck*1e9)*Sys.g;
orbFactor = +bmagn/(planck*1e9)*Sys.orf;
nucFactor = -nmagn/(planck*1e9)*Sys.gn.*Sys.gnscale;

% Loop over all angular momenta (electron spins, nuclear spins, orbitals) selected
for i = idx
  
  if i<=nElectrons
    
    % Electron spin
    if Sys.fullg
      g = elFactor((i-1)*3+(1:3),:);
    else
      g = diag(elFactor(i,:));
    end
    % Transform g matrix to molecular frame
    R_M2g = erot(Sys.gFrame(i,:)); % mol frame -> g frame
    R_g2M = R_M2g.'; % g frame -> mol frame
    g = R_g2M*g*R_g2M.';
    % Build electon Zeeman Hamiltonian in MHz/mT
    for k = 1:3
      Sk = sop(QuantumNumbers,[i,k],'sparse');
      ZxM = ZxM + g(1,k)*Sk;
      ZyM = ZyM + g(2,k)*Sk;
      ZzM = ZzM + g(3,k)*Sk;
    end
    
  elseif i<=nElectrons+nNuclei
    
    % Nuclei, with isotropic gn and chemical shielding (CS) tensor sigma
    iNuc = i-nElectrons;
    if Sys.fullsigma
      sigma = Sys.sigma((iNuc-1)*3+(1:3),:);
    else
      sigma = diag(Sys.sigma(iNuc,:));
    end
    % Transform CS tensor to molecular frame
    ang = Sys.sigmaFrame(iNuc,:);
    if any(ang)
      R_M2CS = erot(ang); % mol frame -> CS frame
      R_CS2M = R_M2CS.'; % CS frame -> mol frame
      sigma = R_CS2M*sigma*R_CS2M.';
    end
    % Build nuclear Zeeman Hamiltonian in MHz/mT
    pre = nucFactor(iNuc);
    for k = 1:3
      Ik = sop(QuantumNumbers,[i,k],'sparse');
      ZxM = ZxM + pre*sigma(1,k)*Ik;
      ZyM = ZyM + pre*sigma(2,k)*Ik;
      ZzM = ZzM + pre*sigma(3,k)*Ik;
    end
    
  else
    
    % Orbital angular momenta, isotropic
    % Build orbital Zeeman Hamiltonian in MHz/mT
    pre = orbFactor(i-nElectrons-nNuclei);
    ZxM = ZxM + pre*sop(QuantumNumbers,[i,1],'sparse');
    ZyM = ZyM + pre*sop(QuantumNumbers,[i,2],'sparse');
    ZzM = ZzM + pre*sop(QuantumNumbers,[i,3],'sparse');
    
  end
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

return
