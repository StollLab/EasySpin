% resfreqs_matrix Compute resonance frequencies for frequency-swept EPR
%
%   ... = resfreqs_matrix(Sys,Exp)
%   ... = resfreqs_matrix(Sys,Exp,Opt)
%   [Pos,Int] = resfreqs_matrix(...)
%   [Pos,Int,Wid] = resfreqs_matrix(...)
%   [Pos,Int,Wid,Trans] = resfreqs_matrix(...)
%
%   Computes frequency-domain EPR line positions, intensities and widths using
%   matrix diagonalization.
%
%   Input:
%    Sys: spin system structure
%    Exp: experimental parameters
%      Field               static field, in mT
%      Range               frequency sweep range, [numin numax], in GHz
%      CenterField         frequency sweep range, [center sweep], in GHz
%      Temperature         temperature, in K
%      SampleFrame         Nx3 array of Euler angles (in radians) for sample/crystal orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%      Mode                excitation mode: 'perpendicular', 'parallel', {k_tilt alpha_pol}
%    Opt: additional computational options
%      Verbosity           level of detail of printing; 0, 1, 2
%      Transitions         nx2 array of level pairs
%      Threshold           cut-off for transition intensity, between 0 and 1
%      Hybrid              0 or 1, switches hybrid mode off or on
%      HybridCoreNuclei    for hybrid mode, nuclei to include in exact core
%      Sites               list of crystal sites to include (default []: all)
%
%   Output:
%    Pos     line positions (in mT)
%    Int     line intensities
%    Wid     Gaussian line widths, full width half maximum (FWHM)
%    Trans   list of transitions included in the computation

function varargout = resfreqs_matrix(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% General
%---------------------------------------------------------------------
% Assert correct Matlab version
warning(chkmlver);

% Guard against wrong number of input or output arguments.
if nargin<1, error('Please supply a spin system as first parameter.'); end
if nargin>3, error('Too many input arguments, the maximum is three.'); end
%Initialize options structure to zero if not given.
if nargin<2, Exp = struct; end
if nargin<3, Opt = struct; end
if isempty(Opt), Opt = struct; end

if nargout<0, error('Not enough output arguments.'); end
if nargout>4, error('Too many output arguments.'); end

if isempty(Exp)
  Exp = struct;
end
if isempty(Opt)
  Opt = struct;
end

if ~isstruct(Sys)
  error('First input argument (Sys) must be a structure!');
end
if ~isstruct(Exp)
  error('Second input argument (Exp) must be a structure!');
end
if ~isstruct(Opt)
  error('Third input argument (Opt) must be a structure!');
end


% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;

% Process Spin system.
%---------------------------------------------------------------------
[Sys,err] = validatespinsys(Sys);
error(err);

DefaultSystem.lw = 0;
DefaultSystem.A=0;
DefaultSystem.HStrain = [0 0 0];
DefaultSystem.gStrain = [0 0 0];
DefaultSystem.AStrain = [0 0 0];
DefaultSystem.DStrain = 0;
DefaultSystem.gAStrainCorr = +1;

Sys = adddefaults(Sys,DefaultSystem);

if (numel(Sys.gAStrainCorr)~=1) || ~isnumeric(Sys.gAStrainCorr) || ...
    (Sys.gAStrainCorr==0) || ~isfinite(Sys.gAStrainCorr)
  error('Sys.gAStrainCorr must be a single number, either +1 or -1.');
end
Sys.gAStrainCorr = sign(Sys.gAStrainCorr);

if Sys.nElectrons>1
  if any(Sys.gStrain(:)) || any(Sys.AStrain(:))
    error('Cannot use D or g/A strain in spin system with more than one electron spin.');
  end
end

if any(Sys.gStrain(:)) || any(Sys.AStrain(:))
  gFull = size(Sys.g,1)==3*numel(Sys.S);
  if gFull
    error('gStrain and AStrain are not allowed when full g matrices are given!');
  end
  if any(Sys.DStrain)
    error('D strain and g/A strain cannot be used at the same time.');
  end
end

if any(Sys.DStrain(:)) && any(Sys.DFrame(:))
  error('D stain cannot be used with tilted D tensors.');
end

if any( strncmp(fieldnames(Sys),'ZB',2))
  higherOrder = 1;
else
  higherOrder = 0;
end

% Process Parameters.
%---------------------------------------------------------------------
DefaultExp.Range = NaN;
DefaultExp.Field = NaN;
DefaultExp.CenterSweep = NaN;
DefaultExp.Temperature = NaN;
DefaultExp.mwMode = 'perpendicular';

DefaultExp.SampleFrame = [0 0 0];
DefaultExp.CrystalSymmetry = 1;
DefaultExp.MolFrame = [0 0 0];
DefaultExp.SampleRotation = [];

Exp = adddefaults(Exp,DefaultExp);

% Check for obsolete fields
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame/Exp.MolFrame instead.');
end
if isfield(Exp,'Mode')
  error('Exp.Mode is no longer supported. Use Exp.mwMode instead.');
end

if isnan(Exp.Field)
  Exp.Field = 0.0;
  logmsg(1,'Exp.Field is missing, assuming 0.0 mT');
end

if ~isnan(Exp.CenterSweep)
  if ~isnan(Exp.Range)
    logmsg(1,'Using Experiment.CenterSweep and ignoring Experiment.Range.');
  end
  Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
  Exp.Range = max(Exp.Range,0);
end

if isnan(Exp.Range), Exp.Range = []; end
if any(diff(Exp.Range)<=0) || any(~isfinite(Exp.Range)) || ~isreal(Exp.Range) || any(Exp.Range<0)
  error('Exp.Range is not valid!');
end

Exp.Range = Exp.Range*1e3; % GHz -> MHz, for comparison with Pdat

% Determine excitation mode
[xi1,xik,nB1_L,nk_L,nB0_L,mwmode] = p_excitationgeometry(Exp.mwMode);

% Photoselection
if ~isfield(Exp,'lightBeam'), Exp.lightBeam = ''; end
if ~isfield(Exp,'lightScatter'), Exp.lightScatter = 0; end

usePhotoSelection = ~isempty(Exp.lightBeam) && Exp.lightScatter<1;

if usePhotoSelection
  if ~isfield(Sys,'tdm') || isempty(Sys.tdm)
    error('To include photoselection weights, Sys.tdm must be given.');
  end
  if ischar(Exp.lightBeam)
    k = [0;1;0]; % beam propagating along yL
    switch Exp.lightBeam
      case 'perpendicular'
        alpha = -pi/2; % gives E-field along xL
      case 'parallel'
        alpha = pi; % gives E-field along zL
      case 'unpolarized'
        alpha = NaN; % unpolarized beam
      otherwise
        error('Unknown string in Exp.lightBeam. Use '''', ''perpendicular'', ''parallel'' or ''unpolarized''.');
    end
    Exp.lightBeam = {k alpha};
  else
    if ~iscell(Exp.lightBeam) || numel(Exp.lightBeam)~=2
      error('Exp.lightBeam should be a 2-element cell {k alpha}.')
    end
  end
end

if ~isfield(Opt,'Sites'), Opt.Sites = []; end

% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt);


% Options parsing and setting.
%---------------------------------------------------------------------

% documented fields
DefaultOptions.Transitions = [];  % list of transitions to include
DefaultOptions.Threshold = 1e-4;  % cutoff threshold for transition pre-selection
DefaultOptions.Hybrid = 0;
DefaultOptions.HybridCoreNuclei = [];

% undocumented fields
DefaultOptions.FuzzLevel = 1e-10;
DefaultOptions.Intensity = 1;
DefaultOptions.Sparse = false;

% Threshold for selecting nuclei for perturbational treatment
DefaultOptions.HybridHFIThreshold = 0.02;
% Threshold for intensities in splitting patterns
DefaultOptions.HybridIntThreshold = 0.005;
% 1 if NQI and NZI of the perturbational nuclei should be removed
DefaultOptions.HybridOnlyHFI = 0;
%:TODO: Remove Options.HybridOnlyHFI
% Remove this option, since it works only in
% strong coupling case. Mims matrix M is (almost) diagonal in
% this case, but also in the limit of weak coupling. Better to
% do full treatment for a general M all the time.

Opt = adddefaults(Opt,DefaultOptions);

if isfield(Opt,'Method')
  Opt.Hybrid = strcmp(Opt.Method,'hybrid');
end

t_ = Opt.Threshold;
if any(~isreal(t_)) || numel(t_)>2 || any(t_<0) || any(t_>=1)
  error('Options.Threshold must be a number >=0 and <1.');
end
preSelectionThreshold = Opt.Threshold(1);
if isscalar(Opt.Threshold)
  postSelectionThreshold = preSelectionThreshold;
else
  postSelectionThreshold = Opt.Threshold(2);
end

StrainsPresent = any([Sys.HStrain(:); Sys.DStrain(:); Sys.gStrain(:); Sys.AStrain(:)]);
computeStrains = StrainsPresent && (nargout>2);
computeIntensities = ((nargout>1) & Opt.Intensity);


% Preparing kernel and perturbing system Hamiltonians.
%-----------------------------------------------------------------------
logmsg(1,'- Preparations');

if Opt.Sparse
  warning('off','MATLAB:eigs:TooManyRequestedEigsForRealSym');
  logmsg(1,'  using sparse matrices');
else
  logmsg(1,'  using full matrices');
end

% Perturbational treatment of SHF nuclei
if Sys.nNuclei>=1 && Opt.Hybrid
  if Opt.Sparse
    error('Cannot use sparse matrices in hybrid mode.');
  end 
  if any(Opt.HybridCoreNuclei>Sys.nNuclei)
    error('At least one entry in Opt.HybridCoreNuclei is incorrect. There are only %d nuclei.',Sys.nNuclei);
  end

  perturbNuclei = true(1,Sys.nNuclei);
  perturbNuclei(Opt.HybridCoreNuclei) = false;  
  idxPerturbNuclei = find(perturbNuclei);
  % :TODO: Allow 1st-order PT only if (2nd-order) error smaller than field increment.
  nPerturbNuclei = numel(idxPerturbNuclei);
  if nPerturbNuclei==0
    str_ = 'none';
  else
    str_ = sprintf('%d ',idxPerturbNuclei);
  end
  logmsg(1,['  nuclei with first-order perturbation treatment: ' str_]);
  
  % Remove perturbational nuclei from core system
  CoreSys = nucspinrmv(Sys,idxPerturbNuclei);
  CoreSys.processed = false;
  CoreSys = rmfield(CoreSys,'lwpp');
  [CoreSys,err] = validatespinsys(CoreSys);
  error(err);
  
  % Prepare terms for nuclear Hamiltonians
  for iiNuc = nPerturbNuclei:-1:1
    iNuc = idxPerturbNuclei(iiNuc);
    I = Sys.I(iNuc);
    [Ix,Iy,Iz] = sop(I,'x','y','z');
    
    % Hyperfine interaction
    for iElectron = 1:Sys.nElectrons
      idxE = 3*(iElectron-1)+(1:3);
      if Sys.fullA
        A = Sys.A(3*(iNuc-1)+(1:3),idxE);
      else
        A = Sys.A(iNuc,idxE);
        R = eye(3);
        if isfield(Sys,'AFrame')
          R = erot(Sys.AFrame(iNuc,idxE)).'; % A frame -> molecular frame
        end
        A = R*diag(A)*R.';
      end
      Hhfi(iElectron,iiNuc).x = A(1,1)*Ix + A(1,2)*Iy + A(1,3)*Iz;
      Hhfi(iElectron,iiNuc).y = A(2,1)*Ix + A(2,2)*Iy + A(2,3)*Iz;
      Hhfi(iElectron,iiNuc).z = A(3,1)*Ix + A(3,2)*Iy + A(3,3)*Iz;
    end
    
    if ~Opt.HybridOnlyHFI
      % Nuclear Zeeman interaction
      prefactor = -nmagn/planck/1e9*Sys.gn(iNuc);
      Hzeem(iiNuc).x = prefactor*Ix;
      Hzeem(iiNuc).y = prefactor*Iy;
      Hzeem(iiNuc).z = prefactor*Iz;
      % Nuclear quadrupole interaction
      Hquad{iiNuc} = 0;
      if I>=1
        Q = [0 0 0];
        R = eye(3);
        if isfield(Sys,'Q'), Q = Sys.Q(iNuc,:); end
        if isfield(Sys,'QFrame')
          R = erot(Sys.QFrame(iNuc,:)).'; % Q frame -> molecular frame
        end
        Q = R*diag(Q)*R.';
        Ivec = {Ix,Iy,Iz};
        for c1 = 1:3
          for c2 = 1:3
            Hquad{iiNuc} = Hquad{iiNuc} + Q(c1,c2)*Ivec{c1}*Ivec{c2};
          end
        end
      end
    end
    
  end
  
  % Components of S vector operator for computing <u|S|u>
  for iEl = Sys.nElectrons:-1:1
    S(iEl).x = sop(CoreSys,[iEl,1]);
    S(iEl).y = sop(CoreSys,[iEl,2]);
    S(iEl).z = sop(CoreSys,[iEl,3]);
  end
  
else
  nPerturbNuclei = 0;
  CoreSys = Sys;
end

% Hamiltonian components for the core system.
if higherOrder
  nCore = hsdim(CoreSys);
  nLevels = nCore;
  nFull = hsdim(Sys);
  nSHFNucStates = nFull/nCore;
else
  if Opt.Sparse
    [kH0,kmuxM,kmuyM,kmuzM] = ham(CoreSys,[],'sparse');
    nLevels = length(kH0);
  else
    [kH0,kmuxM,kmuyM,kmuzM] = ham(CoreSys);
    nLevels = length(kH0);
  end
  nCore = length(kH0);
  nFull = hsdim(Sys);
  nSHFNucStates = nFull/nCore;
end

% Temperature, non-equilibrium populations
computeNonEquiPops = (isfield(Sys,'initState') && ~isempty(Sys.initState));
if computeNonEquiPops

  initState = Sys.initState{1};
  initStateBasis = Sys.initState{2};

  % Check and adapt input dimensions
  nElectronStates = prod(2*Sys.S+1);

  [sz1,sz2] = size(initState);
  if sz1==sz2
    % Density matrix
    if nElectronStates~=sz1 && nCore~=sz1
      error('The density matrix in Sys.initState must have dimensions of nxn with n = %d or %d.',nElectronStates,nCore)
    end
    if numel(initState)==nElectronStates^2 && numel(initState)~=nCore^2
      initState = kron(initState,eye(nCore/nElectronStates))/(nCore/nElectronStates);
    end
  else
    % Vector of populations
    if numel(initState)~=nElectronStates && numel(initState)~=nCore
      error('The population vector in Sys.initState must have %d or %d elements.',nElectronStates,nCore);
    end
    initState = initState(:);
    if numel(initState)==nElectronStates && numel(initState)~=nCore
      initState = kron(initState,ones(nCore/nElectronStates,1))/(nCore/nElectronStates);
    end
  end
  
  computeBoltzmannPopulations = false;
elseif isempty(Exp.Temperature)
  computeBoltzmannPopulations = false;
else
  if numel(Exp.Temperature)~=1
    error('If given, Exp.Temperature must be a single number.');
  end
  if isinf(Exp.Temperature)
    error('If given, Exp.Temperature must be a finite value.');
  end
  computeBoltzmannPopulations = ~isnan(Exp.Temperature);
end

% Add slight numerical noise to non-zero elements in the Hamiltonian to break
% possible degeneracies. Apply if there are more than one electrons or nuclei.
% This is a very crude workaround to prevent numerical issues due to degeneracies.
% It probably adds noise in a lot of situations where it is not necessary.
if Opt.FuzzLevel>0 && ~higherOrder && (CoreSys.nNuclei>1 || CoreSys.nElectrons>1) && ~computeNonEquiPops
  noise = 1 + Opt.FuzzLevel*(2*rand(size(kH0))-1);
  noise = (noise+noise.')/2; % make sure it's Hermitian
  kH0 = kH0.*noise;
  kmuxM = kmuxM.*noise;
  kmuyM = kmuyM.*noise;
  kmuzM = kmuzM.*noise;
end

if nPerturbNuclei>0
  logmsg(1,'  core system with %d spins and %d states',numel(spinvec(CoreSys)),nCore);
  logmsg(1,'  first-order perturbational nuclei with %d states',nSHFNucStates');
else
  if CoreSys.nNuclei>0
    logmsg(1,'  full treatment of all nuclei');
  end
end


% Spin-polarized systems: precompute zero-field energies, states, populations
if computeNonEquiPops && strcmp(initStateBasis,'zerofield')
  
  % Pre-compute zero-field energies and eigenstates
  if higherOrder
    [ZFStates,ZFEnergies] = eig(ham(CoreSys, zeros(1,3)));
  else
    if Opt.Sparse
      [ZFStates,ZFEnergies] = eigs(kH0,length(kH0));
    else
      [ZFStates,ZFEnergies] = eig(kH0);
    end
  end
  [ZFEnergies,idx] = sort(real(diag(ZFEnergies)));
  ZFStates = ZFStates(:,idx);
  % Check for degeneracies and issue error
  if numel(unique(ZFEnergies))~=numel(ZFEnergies)
    error(['Degenerate energy levels detected at zero-field. This prevents unambiguous assignment of ' ...
           'the provided sublevel populations to the zero-field states. Please provide the non-equilibrium' ...
           'state using the full density matrix. See documentation for details.'])
  end

  if isvector(initState)
    % Convert population vector to density matrix
    initState = ZFStates*diag(initState)*ZFStates';
  else
    % Convert density matrix in zero-field basis to uncoupled basis
    initState = ZFStates*initState*ZFStates';
  end

else
  %ZFEnergies = sort(real(eig(kF)));
end


%=======================================================================
% Construction of transition list
%=======================================================================

UserTransitions = ~isempty(Opt.Transitions);
if UserTransitions
  if ischar(Opt.Transitions)
    if strcmp(Opt.Transitions,'all')
      nSStates = prod(2*CoreSys.S+1)*prod(2*CoreSys.L+1);
      logmsg(1,'  using all %d transitions',nSStates*(nSStates-1)/2);
      Transitions = nchoosek(1:nSStates,2);
    else
      error('Options.Transitions must be ''all'' or a nx2 array of enery level indices.');
    end
  else
    % User-specified list of transitions.
    logmsg(1,'  using %d user-specified transitions',size(Opt.Transitions,1));
    % Guarantee that lower index comes first (gives later u < v).
    if size(Opt.Transitions,2)~=2
      error('Options.Transitions must be a nx2 array of energy level indices.');
    end
    Transitions = sort(Opt.Transitions,2);
  end
  
else
  % Automatic compilation: include all level pairs
  Transitions = nchoosek(1:nCore,2);
end

% Terminate if the transition list is empty.
if isempty(Transitions)
  error('No transitions selected! Decrease Opt.Threshold.');
end

if any(Transitions(:)>nCore)
  error('Level index in Options.Transitions is out of range.');
end

% Compute indices and variables used later in the algorithm
u = Transitions(:,1);
v = Transitions(:,2); % u < v
nTransitions = length(u);

% Diagnostic display.
logmsg(1,'  %d transitions pre-selected',nTransitions);

%=======================================================================
% Line width preparations
%=======================================================================
logmsg(1,'- Broadenings');
if computeStrains
  logmsg(1,'  using strains');
  
  % D strain
  %-----------------------------------------------
  [useDStrain,dHdD,dHdE] = getdstrainops(CoreSys);
  
  % g-A strain
  %-------------------------------------------------
  % g strain tensor is taken to be along the g tensor itself.
  usegStrain = any(CoreSys.gStrain(:));
  simplegStrain = CoreSys.nElectrons==1;
  if usegStrain
    logmsg(1,'  g strain present');
    for iEl = CoreSys.nElectrons:-1:1
      gStrainMatrix{iEl} = diag(CoreSys.gStrain(iEl,:)./CoreSys.g(iEl,:));
      if any(CoreSys.gFrame(iEl,:))
        R_g2M = erot(CoreSys.gFrame(iEl,:)).'; % g frame -> molecular frame
        gStrainMatrix{iEl} = R_g2M*gStrainMatrix{iEl}*R_g2M.';
      end
    end
    if ~simplegStrain
      logmsg(1,'  multiple g strains present');
      for iEl = CoreSys.nElectrons:-1:1
        kSxM{iEl} = sop(CoreSys,[iEl,1]);
        kSyM{iEl} = sop(CoreSys,[iEl,2]);
        kSzM{iEl} = sop(CoreSys,[iEl,3]);
      end
    end
  else
    for iEl = CoreSys.nElectrons:-1:1
      gStrainMatrix{iEl} = 0;
    end
  end
  
  useAStrain = (CoreSys.nNuclei>0) && any(CoreSys.AStrain(:));
  if useAStrain
    if isfield(CoreSys,'AFrame')
      R = erot(CoreSys.AFrame(1,:)).'; % A frame -> molecular frame
    else
      R = eye(3);
    end
    
    Ix_ = R(1,1)*sop(CoreSys,[2,1])+R(2,1)*sop(CoreSys,[2,2])+R(3,1)*sop(CoreSys,[2,3]);
    Iy_ = R(1,2)*sop(CoreSys,[2,1])+R(2,2)*sop(CoreSys,[2,2])+R(3,2)*sop(CoreSys,[2,3]);
    Iz_ = R(1,3)*sop(CoreSys,[2,1])+R(2,3)*sop(CoreSys,[2,2])+R(3,3)*sop(CoreSys,[2,3]);
    
    Sx_ = R(1,1)*sop(CoreSys,[1,1])+R(1,2)*sop(CoreSys,[1,2])+R(1,3)*sop(CoreSys,[1,3]);
    Sy_ = R(2,1)*sop(CoreSys,[1,1])+R(2,2)*sop(CoreSys,[1,2])+R(2,3)*sop(CoreSys,[1,3]);
    Sz_ = R(3,1)*sop(CoreSys,[1,1])+R(3,2)*sop(CoreSys,[1,2])+R(3,3)*sop(CoreSys,[1,3]);
    
    dHdAx = CoreSys.AStrain(1)*Ix_*Sx_;
    dHdAy = CoreSys.AStrain(2)*Iy_*Sy_;
    dHdAz = CoreSys.AStrain(3)*Iz_*Sz_;
    
    clear Ix_ Iy_ Iz_ Sx_ Sy_ Sz_
  end
  if any(CoreSys.HStrain), logmsg(2,'  ## using H strain'); end
  if usegStrain, logmsg(2, ' ## using g strain'); end
  if useAStrain, logmsg(2,'  ## using A strain'); end
  if useDStrain, logmsg(2,'  ## using D strain'); end
  
else
  logmsg(1,'  no strains specified',nTransitions);
end


%=======================================================================
%                  DATA GENERATION OVER ORIENTATIONS
%=======================================================================

% Pre-allocations and initializations.
%-----------------------------------------------------------------------
% Pre-allocations are done only when the arrays are needed.

startTime = cputime;
logstr = '';

% Calculate transition rates over all orientations (fixed field!).
Pdat = zeros(nTransitions,nOrientations);

if computeIntensities
  Idat = zeros(nTransitions,nOrientations);
else
  Idat = [];
end

if computeStrains
  Wdat = zeros(nTransitions,nOrientations);
else
  Wdat = [];
end

for iOri = 1:nOrientations
  if EasySpinLogLevel>=1
    if iOri>1
      remainingTime = (cputime-startTime)/(iOri-1)*(nOrientations-iOri+1);
      backspace = repmat(sprintf('\b'),1,numel(logstr));
      hours = fix(remainingTime/3600);
      minutes = fix(remainingTime/60 - 60*hours);
      seconds = remainingTime - 3600*hours - 60*minutes;
      logstr = sprintf('  %d/%d orientations, remaining time %02d:%02d:%0.1f\n', ...
        iOri, nOrientations, hours, minutes, seconds);
      if EasySpinLogLevel==1, fprintf(backspace); end
      fprintf(logstr);
    else
      if nOrientations>1
        logstr = sprintf('  1/%d orientations, remaining time unknown\n',nOrientations);
        fprintf(logstr);
      end
    end
  end
  
  % Set up Hamiltonians for 3 lab principal directions
  %-----------------------------------------------------
  [xLab,yLab,zLab] = erot(Orientations(iOri,:),'rows');
  if higherOrder
    [Vs,E] = gethamdata_hO(Exp.Field,zLab, CoreSys,Opt.Sparse, [], nLevels);
    if Opt.Sparse, sp = 'sparse'; else, sp = ''; end
    [g0e{1},g0e{2},g0e{3}] = ham_ez(CoreSys,[],sp);
    [g0n{1},g0n{2},g0n{3}] = ham_nz(CoreSys,[],sp);
    [g0o{1},g0o{2},g0o{3}] = ham_oz(CoreSys,[],sp);
    for k = 1:3
      g0{k} = g0e{k} + g0n{k} + g0o{k};
    end
    g1 = ham_ezho(CoreSys,[],[],sp,1);
    for n =3:-1:1
      kGM{n} = g0{n}-g1{1}{n};
    end
    % z laboratoy axis: external static field
    kmuzL = zLab(1)*kGM{1} + zLab(2)*kGM{2} + zLab(3)*kGM{3};
    % x laboratory axis: B1 excitation field
    kmuxL = xLab(1)*kGM{1} + xLab(2)*kGM{2} + xLab(3)*kGM{3};
    % y laboratory vector: needed for integration over all B1 field orientations.
    kmuyL = yLab(1)*kGM{1} + yLab(2)*kGM{2} + yLab(3)*kGM{3};
  else
    % z laboratoy axis: external static field
    kmuzL = zLab(1)*kmuxM + zLab(2)*kmuyM + zLab(3)*kmuzM;
    % x laboratory axis: B1 excitation field
    kmuxL = xLab(1)*kmuxM + xLab(2)*kmuyM + xLab(3)*kmuzM;
    % y laboratory vector: needed for integration over all B1 field orientations.
    kmuyL = yLab(1)*kmuxM + yLab(2)*kmuyM + yLab(3)*kmuzM;
    
    [Vs,E] = gethamdata(Exp.Field, kH0, kmuzL, [], nLevels);
  end
  Pdat(:,iOri) = E(v) - E(u);
  
  % Calculate intensities if requested
  if computeIntensities

    % Calculate photoselection weight if needed
    if usePhotoSelection
      k = Exp.lightBeam{1};  % propagation direction
      alpha = Exp.lightBeam{2};  % polarization angle
      if averageOverChi
        ori = Orientations(iOri,1:2);  % omit chi
      else
        ori = Orientations(iOri,1:3);
      end
      photoWeight = photoselect(Sys.tdm,ori,k,alpha);
      % Add isotropic contribution (from scattering)
      photoWeight = (1-Exp.lightScatter)*photoWeight + Exp.lightScatter;
    else
      photoWeight = 1;
    end

    % Compute quantum-mechanical transition rate
    for iTrans = nTransitions:-1:1
      
      U = Vs(:,u(iTrans)); % lower-energy state (u)
      V = Vs(:,v(iTrans)); % higher-energy state (v, Ev>Eu)
      mu_L = [V'*kmuxL*U; V'*kmuyL*U; V'*kmuzL*U]; % magnetic transition dipole moment, in lab frame
      if averageOverChi
        if mwmode.linearpolarizedMode
          TransitionRate = ((1-xi1^2)*norm(mu_L)^2+(3*xi1^2-1)*abs(nB0_L.'*mu_L)^2)/2;
        elseif mwmode.unpolarizedMode
          TransitionRate = ((1+xik^2)*norm(mu_L)^2-(3*xik^2-1)*abs(nB0_L.'*mu_L)^2)/4;
        elseif mwmode.circpolarizedMode
          TransitionRate = ((1+xik^2)*norm(mu_L)^2-(3*xik^2-1)*abs(nB0_L.'*mu_L)^2)/2 - ...
            mwmode.circSense*xik*(nB0_L.'*cross(1i*mu_L,conj(mu_L)));
        end
      else
        if mwmode.linearpolarizedMode
          TransitionRate = abs(nB1_L.'*mu_L)^2;
        elseif mwmode.unpolarizedMode
          TransitionRate = (norm(mu_L)^2-abs(nk_L.'*mu_L)^2)/2;
        elseif mwmode.circpolarizedMode
          TransitionRate = (norm(mu_L)^2-abs(nk_L.'*mu_L)^2) - ...
            mwmode.circSense*(nk_L.'*cross(1i*mu_L,conj(mu_L)));
        end
      end
      if abs(TransitionRate)<1e-12, TransitionRate = 0; end
      TransitionRates(iTrans) = TransitionRate;

    end
    
    TransitionRates = abs(TransitionRates);
    
    
    % Compute polarizations if temperature is given.
    if computeBoltzmannPopulations
      
      Populations = ones(nCore,1);
      
      % Pre-factor for thermal equilibrium populations computations.
      BoltzmannPreFactor = -1e6*planck/boltzm/Exp.Temperature; % MHz^-1
      for iState = 1:nCore
        dE = E(iState) - E(1);
        if dE<1e-10, dE = 0; end % assure we recognize degenerate state even if numerically non-degenerate
        Populations(iState) = exp(BoltzmannPreFactor*dE);
      end
      Populations(isnan(Populations)) = 1;
      Populations = Populations/sum(Populations);
      
      Polarization = Populations(u) - Populations(v);
      if nPerturbNuclei>0
        Polarization = Polarization/prod(2*Sys.I+1);
      end
      
    elseif computeNonEquiPops
      switch initStateBasis
        case 'eigen'
          for iState = 1:nCore
            Populations(iState) = initState(iState,iState);
          end
        otherwise
          for iState = 1:nCore
            Populations(iState) = Vs(:,iState)'*initState*Vs(:,iState);
          end
      end
      Polarization = real(Populations(u) - Populations(v));
      if nPerturbNuclei>0
        Polarization = Polarization/prod(2*Sys.I+1);
      end
       
    else
      % no temperature given
      % same polarization for each electron transition
      %Polarization = Polarization/prod(2*Sys.S+1); % needed to make consistent with high-temp limit
      Polarization = 1/prod(2*Sys.I+1);
    end
    Idat(:,iOri) = TransitionRates(:).*Polarization(:)*photoWeight;
    
  end
  
  % Calculate width if requested.
  %--------------------------------------------------
  if computeStrains
    LineWidthSquared = CoreSys.HStrain.^2*zLab.^2;
    for iTrans = 1:nTransitions
      m = @(Op)Vs(:,v(iTrans))'*Op*Vs(:,v(iTrans)) - Vs(:,u(iTrans))'*Op*Vs(:,u(iTrans));
            
      % H strain: Frequency-domain residual width tensor
      LineWidth2 = LineWidthSquared;
      
      % D strain
      if useDStrain
        for iEl = 1:CoreSys.nElectrons
          LineWidth2 = LineWidth2 + abs(m(dHdD{iEl}))^2;
          LineWidth2 = LineWidth2 + abs(m(dHdE{iEl}))^2;
        end
      end
      
      % A strain
      if useAStrain
        LineWidth2 = LineWidth2 + abs(m(dHdAx))^2;
        LineWidth2 = LineWidth2 + abs(m(dHdAy))^2;
        LineWidth2 = LineWidth2 + abs(m(dHdAz))^2;
      end
      
      % g strain
      if usegStrain
        dg2 = (m(kmuzL)*Exp.Field*zLab.'*gStrainMatrix{1}*zLab)^2;
        LineWidth2 = LineWidth2 + abs(dg2);
      end
      Wdat(iTrans,iOri) = sqrt(LineWidth2);
    end
  end
  
  % First-order approximation for nuclei
  %-------------------------------------------------------
  if nPerturbNuclei>0
    for iTrans = 1 :nTransitions
      U = Vs(:,u(iTrans));
      V = Vs(:,v(iTrans));
      % Calculate <S> for both states involved in the transition
      for iEl = Sys.nElectrons:-1:1
        Su(:,iEl) = [U'*S(iEl).x*U; U'*S(iEl).y*U; U'*S(iEl).z*U];
        Sv(:,iEl) = [V'*S(iEl).x*V; V'*S(iEl).y*V; V'*S(iEl).z*V];
      end
      % Build and diagonalize nuclear sub-Hamiltonians
      for iiNuc = 1:nPerturbNuclei
        Hu = 0;
        Hv = 0;
        % Hyperfine (dependent on S)
        for iEl = 1:Sys.nElectrons
          Hu = Hu + Su(1,iEl)*Hhfi(iEl,iiNuc).x + Su(2,iEl)*Hhfi(iEl,iiNuc).y + Su(3,iEl)*Hhfi(iEl,iiNuc).z;
          Hv = Hv + Sv(1,iEl)*Hhfi(iEl,iiNuc).x + Sv(2,iEl)*Hhfi(iEl,iiNuc).y + Sv(3,iEl)*Hhfi(iEl,iiNuc).z;
        end
        % Nuclear Zeeman and quadrupole (independent of S)
        if ~Opt.HybridOnlyHFI
          Hc = Hquad{iiNuc} + Exp.Field*...
            (zLab(1)*Hzeem(iiNuc).x + zLab(2)*Hzeem(iiNuc).y + zLab(3)*Hzeem(iiNuc).z);
          Hu = Hu + Hc;
          Hv = Hv + Hc;
        end
        % hermitianize (important, otherwise eig returns unsorted values)
        Hu = (Hu+Hu')/2;
        Hv = (Hv+Hv')/2;
        [Vii,dEi] = eig(Hu); dEi = diag(dEi);
        [Vfi,dEf] = eig(Hv); dEf = diag(dEf);
        NucTransitionRates = abs(Vfi'*Vii).^2; % the famous Mims matrix M
        % Compute and store all resonance field shifts and amplitude factors.
        % Intensity thresholds are applied later.
        [vidx,uidx] = find(ones(size(NucTransitionRates)));
        pPdatN{iiNuc}(iTrans,iOri,:) = dEf(vidx(:)) - dEi(uidx(:));
        pIdatN{iiNuc}(iTrans,iOri,:) = NucTransitionRates(:);
      end
    end
  end
  %-------------------------------------------------------

end


logmsg(2,'  ## %2d resonances total from %d level pairs',size(Pdat,1),nTransitions);

% Transitions post-selection
%=======================================================================

% Remove resonances with max intensity below post-selection threshold
%-----------------------------------------------------------------------

idxWeakResonances = [];
if computeIntensities && ~UserTransitions
  maxIntensity = max(abs(Idat),[],2); % maximum over orientations
  absoluteMaxIntensity = max(maxIntensity);
  idxWeakResonances = find(maxIntensity<absoluteMaxIntensity*postSelectionThreshold);
  if ~isempty(idxWeakResonances)
    logmsg(2,'  ## %2d resonances below relative intensity threshold %f',numel(idxWeakResonances),postSelectionThreshold);
  else
    logmsg(2,'  ## all resonances above relative intensity threshold %f',postSelectionThreshold);
  end
else
  logmsg(2,'  ## no intensities computed, no intensity post-selection');
end

if EasySpinLogLevel>=2
  partlyNaN = any(isnan(Pdat),2);
  nChopped = sum(partlyNaN);
  if nChopped>0
    logmsg(2,'  ## %2d transitions partly out of range',nChopped);
  end
  logmsg(2,'  ## %2d transitions fully in range',nTransitions-sum(partlyNaN));
end

% Detect and label out-of-range resonances
%----------------------------------------------
if ~isempty(Exp.Range)
  Pdat(Pdat<Exp.Range(1) | Pdat>Exp.Range(2)) = NaN;
end

idxOutOfRangeResonances = find(all(isnan(Pdat),2));
if ~isempty(idxOutOfRangeResonances)
  logmsg(2,'  ## %2d resonances completely out of range',numel(idxOutOfRangeResonances));
end


% Remove resonances
%----------------------------------------------
idxRmv = [idxOutOfRangeResonances; idxWeakResonances];
if numel(idxRmv)>0
  logmsg(2,'  ## removing %2d resonances (below threshold, out of range)',numel(idxRmv));
  Pdat(idxRmv,:) = [];
  Transitions(idxRmv,:)=[];
  if computeIntensities, Idat(idxRmv,:) = []; end
  if computeStrains, Wdat(idxRmv,:) = []; end
  if nPerturbNuclei>0
    for iiNuc = 1:nPerturbNuclei
      pPdatN{iiNuc}(idxRmv,:,:) = [];
      pIdatN{iiNuc}(idxRmv,:,:) = [];
    end
  end
end
nTransitions = size(Pdat,1);
logmsg(2,'  ## %2d resonances left',nTransitions);


logmsg(1,'  %d significant transitions with resonances in range',nTransitions);

if nTransitions==0
  if ~isempty(Exp.Range)
    fprintf(['WARNING: No resonance frequencies at %g mT between %g and %g GHz.\n'...
      '         Check field value and spectrometer frequency range.\n'],...
      Exp.Field,Exp.Range(1)/1e3,Exp.Range(2)/1e3);
  else
    fprintf('WARNING: No resonances, no frequency range!');
  end
  Output = {[],[],[],[]};
  varargout = Output(1:max(nargout,1));
  return
end

% Assert positive intensities, but only for thermal equilibrium populations
if computeIntensities && (~computeNonEquiPops)
  if any(TransitionRates<0)
    logmsg(-Inf,'*********** Negative intensity encountered in resfields!! Please report! **********');
  end
end
if any(Wdat(:)<0)
  logmsg(-Inf,'*********** Negative width encountered in resfields!! Please report! **************');
end


% Transition post-selection for perturbational nuclei
%---------------------------------------------------------------------
if nPerturbNuclei>0
  logmsg(1,'  transition post-selection on perturbation nuclei');
  
  % Remove low-intensity lines from splittings
  for iiNuc = 1:nPerturbNuclei
    % sum over transitions and orientations
    TotalIntensity = squeeze(sum(sum(pIdatN{iiNuc},1),2));
    idxRmv = TotalIntensity<max(TotalIntensity)*Opt.HybridIntThreshold;
    pPdatN{iiNuc}(:,:,idxRmv) = [];
    pIdatN{iiNuc}(:,:,idxRmv) = [];
    nPertShifts(iiNuc) = size(pPdatN{iiNuc},3);
  end
end

% Compute combined total shifts and amplitudes
%---------------------------------------------------------------------
if nPerturbNuclei>0
  logmsg(1,'  combining shifts from perturbational nuclei');
  
  % Prepare index matrices for line combinations
  for iiNuc = 2:nPerturbNuclei
    [idxComb2{iiNuc},idxComb1{iiNuc}] = ...
      find(ones(nPertShifts(iiNuc),prod(nPertShifts(1:iiNuc-1))));
  end
  
  % Compute total shifts and amplitudes
  for iT = nTransitions:-1:1
    for iO = nOrientations:-1:1
      Pcombined = pPdatN{1}(iT,iO,:);
      Icombined = pIdatN{1}(iT,iO,:);
      for iiNuc = 2:nPerturbNuclei
        Pcombined = Pcombined(idxComb1{iiNuc})  + pPdatN{iiNuc}(iT,iO,idxComb2{iiNuc});
        Icombined = Icombined(idxComb1{iiNuc}) .* pIdatN{iiNuc}(iT,iO,idxComb2{iiNuc});
      end
      pPdat(iT,iO,:) = Pcombined;
      pIdat(iT,iO,:) = Icombined;
    end
  end
  clear pPdatN pIdatN
end

% Combine resonance data with perturbation nuclei
%-------------------------------------------------------------------------------
if nPerturbNuclei>0
  nSubTransitions = size(pPdat,3); % transitions per core system level pair
  expand = @(x) repmat(x,[1,1,nSubTransitions]);
  rearrange = @(x) reshape(permute(x,[1 3 2]),nTransitions*nSubTransitions,[]);
  Pdat = rearrange(expand(Pdat)+pPdat);
  Idat = rearrange(expand(Idat).*pIdat);
  if computeStrains
    if numel(Wdat)>0
      Wdat = rearrange(expand(Wdat));
    end
  end
  Transitions = rearrange(expand(Transitions));
  clear pPdat pIdat
end


% Resonance data summary
%-------------------------------------------------------------------------------
logmsg(2,'  ## resonances min %g MHz, max %g MHz',min(Pdat(:)),max(Pdat(:)));
if computeIntensities
  logmsg(2,'  ## amplitudes min %g, max %g',min(Idat(:)),max(Idat(:)));
end
if computeStrains && numel(Wdat)>0
  logmsg(2,'  ## widths min %g mT, max %g mT',min(Wdat(:)),max(Wdat(:)));
end

% Reshape arrays in the case of crystals with multiple sites
d = dbstack;
pepperCall = numel(d)>1 && strcmp(d(2).name,'pepper');
if ~pepperCall
  if nSites>1
    % Pdat, Idat, Wdat have size [nTransitions, nSites*nOrientations]
    % Resize to [nTransitions*nSites, nOrientations]
    siz = [nTransitions*nSites, numel(Pdat)/nTransitions/nSites];
    Pdat = reshape(Pdat,siz);
    if ~isempty(Idat), Idat = reshape(Idat,siz); end
    if ~isempty(Wdat), Wdat = reshape(Wdat,siz); end
  end
end

% Sort transitions lexicograpically (for each crystal site)
[Transitions, idx] = sortrows(Transitions);
if nSites>1
  idx = idx(:) + (0:nSites-1)*nTransitions;
end
Pdat = Pdat(idx,:);
if ~isempty(Idat), Idat = Idat(idx,:); end
if ~isempty(Wdat), Wdat = Wdat(idx,:); end

% Arrange the output
Output = {Pdat,Idat,Wdat,Transitions};
varargout = Output(1:max(nargout,1));

end
