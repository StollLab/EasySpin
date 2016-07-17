
% zeeman  Zeeman interaction Hamiltonian 
%
%   H = zeeman(SpinSystem, B)
%   H = zeeman(SpinSystem, B, Spins)
%   H = zeeman(SpinSystem, B, Spins, 'sparse')
%   [Zx, Zy, Zz] = zeeman(SpinSystem)
%   [Zx, Zy, Zz] = zeeman(SpinSystem, Spins)
%   [Zx, Zy, Zz] = zeeman(SpinSystem, Spins, 'sparse')
%
%   Returns the Zeeman interaction Hamiltonian for
%   the spins 'Spins' of the spin system 'SpinSystem'.
%
%   Input:
%   - SpinSystem: Spin system structure.
%   - B: Magnetic field vector, in millitesla.
%   - Spins: Vector of spin numbers. For one electron spin: 1
%     is the electron, >=2 are the nuclei. For two electron
%     spins: 1 and 2 electrons, >=3 nuclei, etc. If Spins is
%     omitted, all spins are included.
%   - 'sparse': If given, results returned in sparse format.
%
%   Output:
%   - Zx, Zy, Zz: components of the Zeeman interaction Hamiltonian
%     for the selected spins as defined by Hi=d(H)/d(B_i)
%     i=x,y,z where B_i are the cartesian components of
%     the external field. Units are MHz/mT = 1e9 Hz/T. To get the
%     full Hamiltonian, use H = Zx*B(1)+Zy*B(2)+Zz*B(3), where B is
%     the magnetic field in mT.
%   - H: the Hamiltonian of the Zeeman interaction.

function varargout = zeeman(SpinSystem,varargin)

if (nargin==0), help(mfilename); return; end

if (nargout==2) || (nargout>3), error('Wrong number of output arguments!'); end
singleOutput = nargout<2;

if singleOutput
  if (nargin<1) || (nargin>4), error('Wrong number of input arguments!'); end
  if nargin<2
    error('Field vector (second input, in mT) is missing.')
  else
    Field = varargin{1};
  end
  if (nargin<3), Spins = []; else Spins = varargin{2}; end
  if (nargin<4), opt = ''; else opt = varargin{3}; end
else
  Field = [];
  if (nargin<2), Spins = []; else Spins = varargin{1}; end
  if (nargin<3), opt = ''; else opt = varargin{2}; end
end

if ~ischar(opt)
  error('Last input must be a string, ''sparse''.');
end
sparseResult = strcmp(opt,'sparse');

% Validate spin system
[Sys,err] = validatespinsys(SpinSystem);
error(err);

% Vector of spin quantum numbers
SpinVec = Sys.Spins;

% No 'Spins' specified -> use all
if isempty(Spins), Spins = 1:numel(SpinVec); end

% Validate second argument (Spins)
if any(Spins<1) || any(Spins>length(SpinVec))
  error('Spin indices (2nd input argument) invalid!');
end

% Get number of electrons, nuclei and states
nElectrons = Sys.nElectrons;
nStates = Sys.nStates;

% Initialize Zeeman interaction components to zero
ZxM = sparse(nStates,nStates);
ZyM = sparse(nStates,nStates);
ZzM = sparse(nStates,nStates);

elFactor = bmagn/(planck*1e9)*Sys.g;
nucFactor = -nmagn/(planck*1e9)*Sys.gn;

% Loop over all spins selected
for idx = 1:numel(Spins)
  iSpin = Spins(idx);
  if (iSpin<=nElectrons),  % If it's an electron...
    if Sys.fullg
      g = elFactor((iSpin-1)*3+(1:3),:);
    else
      g = diag(elFactor(iSpin,:));
    end
    % Transform g matrix to molecular frame
    R_M2g = erot(Sys.gFrame(iSpin,:)); % mol frame -> g frame
    R_g2M = R_M2g.'; % g frame -> mol frame
    g = R_g2M*g*R_g2M.';
    % Build EZI Hamiltonian in MHz/mT
    for k = 1:3
      Sk = sop(SpinVec,iSpin,k,'sparse');
      ZxM = ZxM + g(1,k)*Sk;
      ZyM = ZyM + g(2,k)*Sk;
      ZzM = ZzM + g(3,k)*Sk;
    end
  else
    % Nuclei, gn always isotropic
    % Build NZI Hamiltonian in MHz/mT
    pre = nucFactor(iSpin-nElectrons);
    pre = pre * Sys.gnscale(iSpin-nElectrons);
    ZxM = ZxM + pre*sop(SpinVec,iSpin,1,'sparse');
    ZyM = ZyM + pre*sop(SpinVec,iSpin,2,'sparse');
    ZzM = ZzM + pre*sop(SpinVec,iSpin,3,'sparse');
  end
end

if isempty(Field)
  if ~sparseResult
    ZxM = full(ZxM);
    ZyM = full(ZyM);
    ZzM = full(ZzM);
  end
  varargout = {ZxM, ZyM, ZzM};
else
  H = ZxM*Field(1) + ZyM*Field(2) + ZzM*Field(3);
  if ~sparseResult
    H = full(H);
  end
  varargout = {H};
end

return
