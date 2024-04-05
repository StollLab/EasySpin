% resfields  Compute resonance fields for cw EPR 
%
%   ... = resfields(Sys,Exp)
%   ... = resfields(Sys,Exp,Opt)
%   [Pos,Int] = resfields(...)
%   [Pos,Int,Wid] = resfields(...)
%   [Pos,Int,Wid,Trans] = resfields(...)
%
%   Computes cw EPR line positions, intensities and widths using
%   matrix diagonalization and adaptive energy level diagram modelling.
%
%   Input:
%    Sys: spin system structure
%    Exp: experimental parameter settings
%      mwFreq              microwave frequency, in GHz
%      Range               field sweep range, [Bmin Bmax], in mT
%      CenterField         field sweep range, [center sweep], in mT
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

function varargout = resfields(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% General
%---------------------------------------------------------------------
% Assert correct Matlab version
warning(chkmlver);

% Check number of input arguments.
switch nargin
  case 2, Opt = struct;
  case 3
  otherwise
    error('Incorrect number of inputs!');
end

if nargout<0, error('Not enough output arguments.'); end
if nargout>5, error('Too many output arguments.'); end

if isempty(Opt), Opt = struct; end

if ~isstruct(Sys) || ~isstruct(Exp) || ~isstruct(Opt)
  error('SpinSystem, Parameters and Options must be structures!');
end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

% Process Spin system.
%---------------------------------------------------------------------
[Sys,err] = validatespinsys(Sys);
error(err);

DefaultSys.lw = 0;
DefaultSys.HStrain = [0 0 0];
DefaultSys.gStrain = [0 0 0];
DefaultSys.AStrain = [0 0 0];
DefaultSys.DStrain = 0;
DefaultSys.gAStrainCorr = +1;

Sys = adddefaults(Sys,DefaultSys);

if numel(Sys.gAStrainCorr)~=1 || ~isnumeric(Sys.gAStrainCorr) || ...
    Sys.gAStrainCorr==0 || ~isfinite(Sys.gAStrainCorr)
  error('Sys.gAStrainCorr must be a single number, either +1 or -1.');
end
Sys.gAStrainCorr = sign(Sys.gAStrainCorr);

if Sys.nElectrons>1
  if any(Sys.AStrain(:))
    error('AStrain is not supported in spin systems with more than one electron spin.');
  end
end

if any(Sys.gStrain(:)) || any(Sys.AStrain(:))
  gFull = size(Sys.g,1)==3*numel(Sys.S);
  %aFull = size(System.A,1)==3*(1+sum(System.Nucs==','));
  if gFull
    error('gStrain and AStrain are not supported when full g matrices are given!');
  end
  if any(Sys.DStrain)
    error('D strain and g/A strain cannot be used at the same time.');
  end
end

higherOrder = any(strncmp(fieldnames(Sys),'Ham',3));

% Process experimental parameters
%---------------------------------------------------------------------
DefaultExp.mwFreq = NaN;
DefaultExp.Range = NaN;
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

if isnan(Exp.mwFreq), error('Exp.mwFreq is missing!'); end
mwFreq = Exp.mwFreq*1e3; % GHz -> MHz

if ~isnan(Exp.CenterSweep)
  if ~isnan(Exp.Range)
    %logmsg(0,'Using Exp.CenterSweep and ignoring Exp.Range.');
  end
  Exp.Range = Exp.CenterSweep(1) + [-1 1]*Exp.CenterSweep(2)/2;
  if Exp.Range(1)<0
    error('Lower field limit from Exp.CenterSweep cannt be negative.');
  end
end

if isnan(Exp.Range), error('Exp.Range/Exp.CenterSweep is missing!'); end
if any(diff(Exp.Range)<=0) || any(~isfinite(Exp.Range)) || ~isreal(Exp.Range)
  error('Exp.Range is not valid!');
end
if any(Exp.Range<0)
  error('Negative magnetic fields in Exp.Range are not possible.');
end

% Determine excitation mode
[xi1,xik,nB1_L,nk_L,nB0_L,mwmode] = p_excitationgeometry(Exp.mwMode);

if ~isfield(Opt,'Sites'), Opt.Sites = []; end

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

% Process crystal orientations, crystal symmetry, and frame transforms
[angles_M2L,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt);

% Options parsing and setting.
%---------------------------------------------------------------------
% obsolete fields
obsoleteOptions = {''};
for iOpt = 1:numel(obsoleteOptions)
  if isfield(Opt,obsoleteOptions{iOpt})
    error('Options.%s is obsolete. Please remove from code!',obsoleteOptions{iOpt});
  end
end

if isfield(Opt,'Perturb')
  error('Options.Perturb is obsolete. Use Opt.Method=''perturb'' or Opt.Method=''hybrid'' instead.');
end

% Fields with changed format: formerly strings/documented, now {0,1}/undocumented
changedFields = {'Intensity','Gradient'};
for iFld = 1:numel(changedFields)
  if isfield(Opt,changedFields{iFld})
    if ischar(Opt.(changedFields{iFld}))
      error('Options.%s is obsolete. Please remove from code!',changedFields{iFld});
    end
  end
end

% documented fields
DefaultOptions.Transitions = [];  % list of transitions to include
DefaultOptions.Threshold = 1e-4;  % cutoff threshold for transition pre-selection
DefaultOptions.Hybrid = 0;
DefaultOptions.HybridCoreNuclei = [];

% undocumented fields
DefaultOptions.TPSGridSize = 4;      % grid size for transition pre-selection
DefaultOptions.TPSGridSymm = 'D2h';  % grid symmetry for transition pre-selection
DefaultOptions.FuzzLevel = 1e-10;
DefaultOptions.Freq2Field = true;
DefaultOptions.maxSegments = 2000;
DefaultOptions.ModellingAccuracy = 2e-6;
DefaultOptions.RediagLimit = 0.95;
DefaultOptions.Sparse = false;
DefaultOptions.nLevels = [];
DefaultOptions.Gradient = 1;
DefaultOptions.Intensity = 1;

% Threshold for selecting nuclei for perturbational treatment
DefaultOptions.HybridHFIThreshold = 0.02;
% Threshold for intensities in splitting patterns
DefaultOptions.HybridIntThreshold = 0.005;
% 1 if NQI and NZI of the perturbational nuclei should be removed
DefaultOptions.HybridOnlyHFI = 0;
%:TODO: Remove Options.HybridOnlyHFI
% Remove this option, doesn't make sense. Works only in
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

IntensitySwitch = Opt.Intensity;
GradientSwitch = Opt.Gradient;

if Opt.Freq2Field~=1 && Opt.Freq2Field~=0
  error('Options.Freq2Field incorrect!');
end
computeFreq2Field = Opt.Freq2Field;

StrainsPresent = any([Sys.HStrain(:); Sys.DStrain(:); Sys.gStrain(:); Sys.AStrain(:)]);
computeStrains = StrainsPresent && (nargout>2);

computeGradient = (computeStrains || (nargout>4)) && GradientSwitch;
computeIntensities = (nargout>1 && IntensitySwitch) || computeGradient;
computeEigenPairs = computeIntensities || computeGradient || computeStrains;

% Preparing kernel and perturbing system Hamiltonians.
%-----------------------------------------------------------------------
logmsg(1,'- Preparations');

if Opt.Sparse
  warning('off','MATLAB:eigs:TooManyRequestedEigsForRealSym');
  logmsg(1,'  using sparse matrices');
else
  logmsg(1,'  using full matrices');
end

CoreSys = Sys;

% HFI splitting at zero field relative to mw frequency
HFIStrength = 0;
if Sys.nNuclei>0
  if Sys.fullA
    for iNuc = Sys.nNuclei:-1:1
      maxHF(iNuc) = max(max(Sys.A((iNuc-1)*3+(1:3))));
    end
  else
    maxHF = max(abs(Sys.A),[],2);
  end
  HFIStrength = maxHF(:).'.*(1/2+nucspin(Sys.Nucs))/mwFreq;
end

% Perturbational treatment of SHF nuclei
if CoreSys.nNuclei>=1 && Opt.Hybrid

  %if iscell(CoreSys.Nucs), CoreSys.Nucs = CoreSys.Nucs{1}; end
  
  if any(Opt.HybridCoreNuclei>CoreSys.nNuclei)
    error('Opt.HybridCoreNuclei is incorrect!');
  end
  perturbNuclei = ones(1,CoreSys.nNuclei);
  perturbNuclei(Opt.HybridCoreNuclei) = 0;
  
  idx = find(perturbNuclei);
  %idx = idx & (HFIStrength<Opt.HybridHFIThreshold);
  % :TODO: Allow 1st-order PT only if (2nd-order) error smaller than field increment.  
  nPerturbNuclei = numel(idx);
  str1 = sprintf('%d ',idx);
  if isempty(str1), str1 = 'none'; end
  logmsg(1,['  nuclei with first-order treatment: ' str1]);

  % Remove perturbational nuclei from core system
  CoreSys = nucspinrmv(CoreSys,idx);
  CoreSys.processed = 0;
  CoreSys = rmfield(CoreSys,'lwpp');
  [CoreSys,err] = validatespinsys(CoreSys);
  error(err);
  
  % Prepare terms for nuclear Hamiltonians
  for iiNuc = nPerturbNuclei:-1:1
    iNuc = idx(iiNuc);
    I = Sys.I(iNuc);
    [Ix,Iy,Iz] = sop(I,'x','y','z');
    nPerturbTransitions(iiNuc) = (2*I+1)^2;
    
    % Hyperfine interaction
    for iElectron = 1:Sys.nElectrons
      idxE = 3*(iElectron-1)+(1:3);
      if Sys.fullA
        A = Sys.A(3*(iNuc-1)+(1:3),idxE);
      else
        A = Sys.A(iNuc,idxE);
        R_A2M = eye(3);
        if isfield(Sys,'AFrame')
          R_A2M = erot(Sys.AFrame(iNuc,idxE)).'; % A frame -> molecular frame
        end
        A = R_A2M*diag(A)*R_A2M.';
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
        R_Q2M = eye(3);
        if isfield(Sys,'Q')
          Q = Sys.Q(iNuc,:);
        end
        if isfield(Sys,'QFrame')
          R_Q2M = erot(Sys.QFrame(iNuc,:)).'; % Q frame -> molecular frame
        end
        Q = R_Q2M*diag(Q)*R_Q2M.';
        Ivec = {Ix,Iy,Iz};
        for c1 = 1:3
          for c2 = 1:3
            Hquad{iiNuc} = Hquad{iiNuc} + Q(c1,c2)*Ivec{c1}*Ivec{c2};
          end
        end
      end
    end

  end
  
  % Components of S vectors for computing <u|S|u>
  for iEl = Sys.nElectrons:-1:1
    S(iEl).x = sop(CoreSys,[iEl,1]);
    S(iEl).y = sop(CoreSys,[iEl,2]);
    S(iEl).z = sop(CoreSys,[iEl,3]);
  end

else
  nPerturbNuclei = 0;
end

% Hamiltonian components for the core system.
if higherOrder
  nCore = hsdim(CoreSys);
  nLevels = hsdim(CoreSys);
  % spin Hamiltonian is calculated later
else
  if Opt.Sparse
    [kH0,kmuxM,kmuyM,kmuzM] = ham(CoreSys,[],'sparse');
    nLevels = length(kH0);
  else
    [kH0,kmuxM,kmuyM,kmuzM] = ham(CoreSys);
    nLevels = length(kH0);
  end
  nCore = length(kH0);
end
nFull = hsdim(Sys);
nSHFNucStates = nFull/nCore;

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
    if sz1~=nElectronStates && sz1~=nCore
      error('The density matrix in Sys.initState must have dimensions of nxn with n = %d or %d.',nElectronStates,nCore)
    end
    if numel(initState)==nElectronStates^2 && nCore>nElectronStates
      initState = kron(initState,eye(nCore/nElectronStates))/(nCore/nElectronStates);
    end
  elseif isvector(initState)
    % Vector of populations
    if numel(initState)~=nElectronStates && numel(initState)~=nCore
      error('The population vector in Sys.initState must have %d or %d elements.',nElectronStates,nCore);
    end
    initState = initState(:);
    if numel(initState)==nElectronStates && numel(initState)~=nCore
      initState = kron(initState,ones(nCore/nElectronStates,1))/(nCore/nElectronStates);
    end
    initState = diag(initState);
  else
    error('Sys.initState must contain either a population vector or a density matrix.');
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
    [ZFStates,ZFEnergies] =  eig(ham(CoreSys,zeros(1,3)));
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
           'the provided sublevel populations to the zero-field states. Please provide the non-equilibrium ' ...
           'state using the full density matrix. See documentation for details.'])
  end
  
  % Convert density matrix in zero-field basis to uncoupled basis
  initState = ZFStates*initState*ZFStates';
  
else
  if higherOrder
    ZFEnergies = eig(ham(CoreSys,zeros(1,3)));
    ZFEnergies = sort(real(ZFEnergies));
  else
    if issparse(kH0)
      ZFEnergies(1) = eigs(kH0,1,-2*max(abs(kH0(:))));
      ZFEnergies(2) = eigs(kH0,1,+2*max(abs(kH0(:))));
    else
      ZFEnergies = sort(real(eig(kH0)));
    end
  end
end

% Check whether looping transitions are possible.
maxZeroFieldSplit = ZFEnergies(end) - ZFEnergies(1);
if maxZeroFieldSplit>1e3
  logmsg(2,'  ## maximum splitting at zero field: %g GHz',maxZeroFieldSplit/1e3);
else
  logmsg(2,'  ## maximum splitting at zero field: %g MHz',maxZeroFieldSplit);
end
LoopFields = maxZeroFieldSplit/mwFreq > 1 - 1e-6;
if LoopFields
  msg = '  looping transitions possible';
else
  msg = '  no looping transitions possible';
end
logmsg(1,msg);


%===============================================================================
% Transition pre-selection
%===============================================================================
% Compose a list of transitions, ie level pairs, for which peak data are
% to be computed. The list is either taken from a user specified list
% in Opt.Transitions, or is determined by an automatic procedure which
% selects the most intense transitions, their number being determined
% by a threshold for the relative transitions rate. The relative transition
% rate for the most intense transition is 1.
%
% The transition rates are determined at center field. As a consequence, this
% pre-selection fails when states vary significantly over the field range, e.g.
% for systems with large zero-field splittings and wide field ranges starting
% at very low fields.

logmsg(1,'- Transition pre-selection');

UserTransitions = ~isempty(Opt.Transitions);
if UserTransitions

  if ischar(Opt.Transitions)
    if strcmp(Opt.Transitions,'all')
      if isempty(Opt.nLevels)
        nStates_ = prod(2*CoreSys.S+1)*prod(2*CoreSys.I+1);
        if isfield(CoreSys,'L')
          nStates_ = nStates_*prod(2*CoreSys.L+1);
        end
      else
        nStates_ = Opt.nLevels;
      end
      logmsg(1,'  using all %d transitions',nStates_*(nStates_-1)/2);
      [u,v] = find(triu(ones(nStates_),1));
      Transitions = sortrows([u,v]);
      postSelectionThreshold = 0;
    else
      error('Options.Transitions must be ''all'' or a nx2 array of enery level indices.');
    end
  else
    % User-specified list of transitions
    logmsg(1,'  using %d user-specified transitions',size(Opt.Transitions,1));
    % Guarantee that lower index comes first (gives later u < v).
    if size(Opt.Transitions,2)~=2
      error('Options.Transitions must be a nx2 array of energy level indices.');
    end
    Transitions = sort(Opt.Transitions,2);
    rmv = any(Transitions>nLevels,2);
    Transitions(rmv,:) = [];
  end

else % Automatic transition pre-selection
  
  nElStates_ = prod(2*CoreSys.S+1)*prod(2*CoreSys.L+1);
  if preSelectionThreshold==0
    logmsg(1,'  selection threshold is zero -> using all transitions');
    TransitionRates = ones(nCore);
  elseif CoreSys.nElectrons==1 && CoreSys.nNuclei==0 ...
      && (~any(CoreSys.L))
    logmsg(1,'  one electron spin and no nuclei -> using all transitions');
    TransitionRates = ones(nCore);
  else
    if nOrientations>1 % if powder or multiple orientations
      % Set a coarse grid, independent of the Hamiltonian symmetry
      logmsg(1,'  selection threshold %g',preSelectionThreshold);
      logmsg(2,'  ## (selection threshold %g, grid size %d, grid symmetry %s)',...
        preSelectionThreshold,Opt.TPSGridSize,Opt.TPSGridSymm);
      TPSgrid = sphgrid(Opt.TPSGridSymm,Opt.TPSGridSize);
      phi = TPSgrid.phi;
      theta = TPSgrid.theta;
      TPSweights = TPSgrid.weights;
    else % single orientation
      phi = angles_M2L(1);
      theta = angles_M2L(2);
      TPSweights = 1;
    end
    
    % Prepare detection operators
    if higherOrder
      if Opt.Sparse
        sp = 'sparse';
        g1 = ham_ezho(CoreSys,[],sp,1);
      else
        sp = '';
        g1 = ham_ezho(CoreSys,[],[],sp,1);
      end
      [g0{1},g0{2},g0{3}] = ham_ez(CoreSys,[],sp);
      if Sys.nNuclei>0
        [mu0n{1},mu0n{2},mu0n{3}] = ham_nz(CoreSys,[],sp);
        for k = 1:3
          g0{k} = g0{k} - mu0n{k};
        end
      end
      ExM = g1{1}{1} + g0{1};
      EyM = g1{1}{2} + g0{2};
      EzM = g1{1}{3} + g0{3};
    else
      ExM = -kmuxM;
      EyM = -kmuyM;
      EzM = -kmuzM;
    end
    
    % Pre-compute trigonometric functions
    st = sin(theta);
    ct = cos(theta);
    sp = sin(phi);
    cp = cos(phi);
    centerB = mean(Exp.Range); % take field at center of scan range

    % Calculate transition rates over all orientations (at fixed field)
    TransitionRates = zeros(nCore);  % preallocate
    for iOri = 1:numel(theta)
      % Determine eigenvectors
      if higherOrder
        [Vecs,~] = gethamdata_hO(centerB,[st(iOri)/sqrt(2)*[1,1],ct(iOri)],CoreSys,Opt.Sparse,[],nLevels);
      else
        kmuzL = st(iOri)*(cp(iOri)*kmuxM + sp(iOri)*kmuyM) + ct(iOri)*kmuzM;
        % Solve eigenproblem
        if Opt.Sparse
          [Vecs,E] = eigs(kH0 - centerB*kmuzL,nCore);
          [~,idx_] = sort(diag(E));
          Vecs = Vecs(:,idx_);
        else
          [Vecs,~] = eig(kH0 - centerB*kmuzL);
        end
      end
      % Calculate transition rate matrix and take the maximum
      ExyM = cp(iOri)*ExM + sp(iOri)*EyM;
      if mwmode.parallelMode
        EzL = st(iOri)*ExyM + ct(iOri)*EzM;
        TransitionRates = max(TransitionRates,TPSweights(iOri) * abs(Vecs'*EzL*Vecs).^2);
      else % perpendicular
        EyL = -sp(iOri)*ExM  + cp(iOri)*EyM;
        ExL =  ct(iOri)*ExyM - st(iOri)*EzM;
        TransitionRates = max(TransitionRates,TPSweights(iOri) * (abs(Vecs'*ExL*Vecs).^2 + abs(Vecs'*EyL*Vecs).^2)/2);
      end
    end
    % Free unused memory
    clear Vecs E idx kmuzL ExM EyM EzM ExyM ExL EyL EzL
  end
  
  % Remove lower triangular part
  idxLowerTriangle = logical(tril(ones(nCore)));
  % Remove nuclear transitions
  if max(HFIStrength)<0.5 && preSelectionThreshold>0
    idxNuclearTransitions = logical(kron(eye(nElStates_),ones(nCore/nElStates_)));
    keepidx = ~idxLowerTriangle & ~idxNuclearTransitions;
  else
    keepidx = ~idxLowerTriangle;
  end
  TransitionRates(~keepidx) = [];
  [u,v] = find(keepidx); % get level indices for transition pairs
  Transitions = [u,v];
  clear keepidx u v idxLowerTriangle idxNuclearTransitions

  % Use threshold for number determination
  nTransitions = sum(TransitionRates>preSelectionThreshold*max(TransitionRates));
  
  % Sort transition rates in descending order
  [~,idx] = sort(TransitionRates,'descend');
  % Select most intense transitions
  Transitions = Transitions(idx(1:nTransitions),:);
  clear TransitionRates idx
  
end

% Terminate if there is a problem with the transition list
if isempty(Transitions)
  error('No transitions selected! Decrease Opt.Threshold.');
end
if any(Transitions(:)>nCore)
  error('Level index in Options.Transitions is out of range.');
end

% Compute indices and variables used later in the algorithm
u = Transitions(:,1);
v = Transitions(:,2); % v > u
nTransitions = length(u);
upTRidx = u + (v-1)*nCore; % Indices into UPPER triangle
Trans = upTRidx; % One-number transition indices

% Diagnostic display
logmsg(1,'  %d transitions pre-selected',nTransitions);

% Now, if the transitions were selected automatically, they are in order of
% descending expeced intensity. If user-specified, their order is unchanged.
%===============================================================================


%=======================================================================
% Line width preparations
%=======================================================================
logmsg(1,'- Broadenings');
simplegStrain = true;
usegStrain = false;
useAStrain = false;
if computeStrains
  logmsg(1,'  using strains');
  
  % Frequency-domain residual width tensor
  %-----------------------------------------------
  HStrain2 = CoreSys.HStrain.^2;
  
  % D strain
  %-----------------------------------------------
  [useDStrain,dHdD,dHdE] = getdstrainops(CoreSys);
  
  % g-A strain
  %-------------------------------------------------
  % g strain tensor is taken to be aligned with the g tensor
  % A strain tensor is taken to be aligned with the A tensor
  % g strain can be specified for each electron spin
  % A strain is limited to the first electron and first nuclear spin
  usegStrain = any(CoreSys.gStrain(:));
  if usegStrain
    logmsg(1,'  g strain present');
    simplegStrain = CoreSys.nElectrons==1;
    for iEl = CoreSys.nElectrons:-1:1
      gStrainMatrix{iEl} = diag(CoreSys.gStrain(iEl,:)./CoreSys.g(iEl,:))*mwFreq; % MHz
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
    for e = CoreSys.nElectrons:-1:1
      gStrainMatrix{e} = zeros(3);
    end
  end
  
  useAStrain = (CoreSys.nNuclei>0) && any(CoreSys.AStrain);
  if useAStrain
    % Transform A strain matrix to molecular frame.
    AStrainMatrix = diag(CoreSys.AStrain);
    if isfield(CoreSys,'AFrame')
      R_A2M = erot(CoreSys.AFrame(1,:)).'; % A frame -> molecular frame
      AStrainMatrix = R_A2M*AStrainMatrix*R_A2M.';
    end
    % Diagonalize Hamiltonian at center field.
    centerB = mean(Exp.Range);
    [Vecs,E] = eig(kH0 - centerB*kmuzM);
    [~,idx] = sort(real(diag(E)));
    Vecs = Vecs(:,idx);
    % Calculate effective mI of nucleus 1 for all eigenstates.
    mI = real(diag(Vecs'*sop(CoreSys,[2,3])*Vecs));
    mITr = mean(mI(Transitions),2);
    % compute A strain array
    AStrainMatrix = reshape(mITr(:,ones(1,9)).',[3,3,nTransitions]).*...
      repmat(AStrainMatrix,[1,1,nTransitions]);
    corr = Sys.gAStrainCorr;
    for e = Sys.nElectrons:-1:1
      gAslw2{e} = (repmat(gStrainMatrix{e},[1,1,nTransitions])+corr*AStrainMatrix).^2;
    end
    clear AStrainMatrix Vecs E idx mI mITr
  else
    for e = Sys.nElectrons:-1:1
      gAslw2{e} = repmat(gStrainMatrix{e}.^2,[1,1,nTransitions]);
    end
  end
  clear gslw
  % gAslw2 = a (cell array of) 3D array with 3x3 strain line-width matrices
  % for each transition piled up along the third dimension.
  
  if any(HStrain2), logmsg(2,'  ## using H strain'); end
  if usegStrain || useAStrain, logmsg(2,'  ## using g/A strain'); end
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

msg = 'positions';
Pdat = NaN(nTransitions,nOrientations);
if nPerturbNuclei>0
  for iiNuc = nPerturbNuclei:-1:1
    pPdatN{iiNuc} = zeros(nTransitions,nOrientations,nPerturbTransitions(iiNuc));
  end
end

if computeIntensities
  Idat = NaN(nTransitions,nOrientations);
  if nPerturbNuclei>0
    for iiNuc = nPerturbNuclei:-1:1
      pIdatN{iiNuc} = zeros(nTransitions,nOrientations,nPerturbTransitions(iiNuc));
    end
  end
  if computeBoltzmannPopulations
    % Pre-factor for thermal equilibrium populations computations
    BoltzmannPreFactor = 1e6*planck/boltzm/Exp.Temperature; % MHz^-1
  end
  msg = [msg ', intensities'];
else
  Idat = [];
end

if computeStrains
  msg = [msg ', strain widths'];
  Wdat = NaN(nTransitions,nOrientations);
else
  Wdat = [];
end

% Pre-allocate Gradient if it is requested.
if computeGradient
  msg = [msg ', gradients'];
  Gdat = zeros(nTransitions,nOrientations);
else
  Gdat = [];
end
logmsg(1,'- Resonance data: computing %s',msg);

idxTr = [(1:nTransitions).'; Inf];

% Preparation for the adaptive iterative bisection
%---------------------------------------------------------
levelAccuracy = mwFreq*Opt.ModellingAccuracy;

M = [2 -2 1 1; -3 3 -2 -1; 0 0 1 0; 1 0 0 0];
ZeroRow = NaN(1,nOrientations);
if higherOrder
  maxSlope = 0;
  for iOri = 1:nOrientations
    [~,~,zLab_M] = erot(angles_M2L(iOri,:),'rows');
    [~,~,der]= gethamdata_hO(Exp.Range(2),zLab_M,CoreSys,Opt.Sparse,[],nLevels);
    maxSlope = max([maxSlope,max(der)]);
  end
else
  if Opt.Sparse
    maxSlope = max(max([eigs(-kmuxM,1) eigs(-kmuyM,1) eigs(-kmuzM,1)]));
    maxSlope = abs(maxSlope);
  else
    maxSlope = max(max([eig(-kmuxM) eig(-kmuyM) eig(-kmuzM)]));
  end
end
nDiagonalizations = 0; % or 4? (1 for F and 3 for maxSlope)
nRediags = 0;
minStateStability = 2;
nMaxSegmentsReached = 0;

logmsg(2,'  ## modelling accuracy %0.3f MHz, max slope %g MHz/mT',levelAccuracy,maxSlope);

% Loop over all given orientations
%-----------------------------------------------------------------------

startTime = cputime;
logstr = '';
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
  % xLab, yLab, zLab represented in the molecular frame M
  [xLab_M,yLab_M,zLab_M] = erot(angles_M2L(iOri,:),'rows');
  
  if ~higherOrder
    % zLab axis: external static field
    kmuzL = zLab_M(1)*kmuxM + zLab_M(2)*kmuyM + zLab_M(3)*kmuzM;
    % xLab axis: mw excitation field
    kmuxL = xLab_M(1)*kmuxM + xLab_M(2)*kmuyM + xLab_M(3)*kmuzM;
    % yLab axis: needed for gradient calculation
    % and the integration over all mw field orientations
    kmuyL = yLab_M(1)*kmuxM + yLab_M(2)*kmuyM + yLab_M(3)*kmuzM;
    if usegStrain && ~simplegStrain
      for e = Sys.nElectrons:-1:1
        kSzL{e} = zLab_M(1)*kSxM{e} + zLab_M(2)*kSyM{e} + zLab_M(3)*kSzM{e};
      end
    end
  end
  if computeStrains
    LineWidthSquared = HStrain2*zLab_M.^2;
  end
  
  % Pre-calculate photoselection weight if needed
  if usePhotoSelection
    k = Exp.lightBeam{1};  % propagation direction
    alpha = Exp.lightBeam{2};  % polarization angle
    if averageOverChi
      ori = angles_M2L(iOri,1:2);  % omit chi
    else
      ori = angles_M2L(iOri,1:3);
    end
    photoWeight = photoselect(Sys.tdm,ori,k,alpha);
    % Add isotropic contribution (from scattering)
    photoWeight = (1-Exp.lightScatter)*photoWeight + Exp.lightScatter;
  else
    photoWeight = 1;
  end
  
  %=============================================================================
  % Build cubic spline model of energy level diagram using iterative bisection
  %-----------------------------------------------------------------------------
  
  nSegments = 1;  % number of segments to start width
  Bknots = linspace(Exp.Range(1),Exp.Range(2),nSegments+1);  % segments with equal width
  
  % Preallocations
  vectors = cell(1,nSegments+1);  % eigenvectors
  E = cell(1,nSegments+1);  % energies
  dEdB = cell(1,nSegments+1);  % dE/dB
  deltaE = cell(1,nSegments+1);  % transition energies for listed transitions
  
  % Calculate energies etc. at all segment boundaries (knots)
  for s = 1:nSegments+1
    if higherOrder
      [vectors{s},E{s},dEdB{s},deltaE{s}] = gethamdata_hO(Bknots(s),zLab_M,CoreSys,Opt.Sparse,Trans,nLevels);
    else
      [vectors{s},E{s},dEdB{s},deltaE{s}] = gethamdata(Bknots(s),kH0,kmuzL,Trans,nLevels);
    end
    nDiagonalizations = nDiagonalizations + 1;
  end
  
  % Iterative bisection to model energy level diagram to desired accuracy
  converged = false(1,nSegments);
  while ~all(converged) && nSegments<Opt.maxSegments
    
    s = find(~converged,1);  % find a segment that's not yet converged
    dB = Bknots(s+1) - Bknots(s);
    if E{s+1}(end)-E{s+1}(1) > mwFreq
      if LoopFields
        transitionPossible = abs((deltaE{s}+deltaE{s+1})/2-mwFreq) <= maxSlope*dB;
      else
        transitionPossible = (deltaE{s}-mwFreq).*(deltaE{s+1}-mwFreq) <= 0;
      end
    else
      transitionPossible = false;
    end
    
    % If transition within segment is not possible, mark segment as converged and move on
    if ~any(transitionPossible)
      converged(s) = true;
      continue
    end
    
    % If transition is possible, diagonalize at midpoint and compute error
    newB = (Bknots(s)+Bknots(s+1))/2;
    if higherOrder
      [Ve,En,dEdBn,dEn] = gethamdata_hO(newB,zLab_M,CoreSys,Opt.Sparse,Trans,nLevels);
    else
      [Ve,En,dEdBn,dEn] = gethamdata(newB,kH0,kmuzL,Trans,nLevels);
    end
    nDiagonalizations = nDiagonalizations + 1;

    % Convergence check - exclude levels that are not involved in resonances
    include = false(1,nCore);
    include([u(transitionPossible) v(transitionPossible)]) = true;
    deviation = 2*((E{s}+E{s+1})/2 + dB/8*(dEdB{s}-dEdB{s+1}) - En);
    deviation = abs(deviation(include));
    converged_ = max(deviation) <= levelAccuracy;
    
    % Bisect segment by adding midpoint
    Bknots = [Bknots(1:s) newB Bknots(s+1:end)];
    E = [E(1:s) {En} E(s+1:end)];
    vectors = [vectors(1:s) {Ve} vectors(s+1:end)];
    dEdB = [dEdB(1:s) {dEdBn} dEdB(s+1:end)];
    deltaE = [deltaE(1:s) {dEn} deltaE(s+1:end)];
    nSegments = nSegments + 1;

    % Mark the two new segments converged or not depending on error
    converged = [converged(1:s-1) converged_ converged_ converged(s+1:end)];
    
  end
  
  if nSegments>=Opt.maxSegments
    nMaxSegmentsReached = nMaxSegmentsReached + 1;
    logmsg(2,'   segmentation terminated, %d segments, max number of segment reached',nSegments);
  else
    logmsg(2,'   segmentation converged, %d segments',nSegments);
  end
  
  
  % Compute eigenvector cross products to determine how strongly eigenvectors
  % change over a segment.
  for s = nSegments:-1:1
    %StateStability(:,s) = abs(diag(Vectors{s+1}'*Vectors{s}));  % slower
    StateStability(:,s) = abs(sum(conj(vectors{s+1}).*vectors{s})).';
  end

  % Cubic polynomial coefficients of the entire spline model of the energy level diagram
  dB = diff(Bknots);
  SplineModelCoeffs = cell(1,nSegments);
  for s = 1:nSegments
    SplineModelCoeffs{s} = M*[E{s}; E{s+1}; dB(s)*dEdB{s}; dB(s)*dEdB{s+1}];
  end
  
  iiTrans = 1;
  for iTrans = 1:nTransitions % run over all level pairs

    % Find first position of transition iTrans in idxTr.
    while idxTr(iiTrans)<iTrans, iiTrans=iiTrans+1; end

    % Loop over all field segments and compute resonance fields
    %-----------------------------------------------------------
    for s = 1:nSegments

      % Construct cubic polynomial coeffs of resonance function Ev-Eu-freq
      coeffs = SplineModelCoeffs{s};
      Evpoly = coeffs(:,v(iTrans));
      Eupoly = coeffs(:,u(iTrans));
      dEpoly = Evpoly - Eupoly;
      dEpoly(4) = dEpoly(4) - mwFreq;

      % Find zeros, first and second derivatives of E(v)-E(u)-mwFreq
      %[Zeros,Diff1,Diff2] = cubicsolve(C,LoopFields);

      % Problem here: cubicsolve should only be called if a resonance
      % is possible. This can be determined as above using maxSlope etc.
      % cubicsolve takes too long to find this out (esp. for LoopFields=1).
      cubicZeros = cubicsolve(dEpoly,LoopFields);
      if isempty(cubicZeros), continue, end

      ResonanceFields = Bknots(s) + dB(s)*cubicZeros;

      % loop over all resonances for transition iTrans in the current segment
      for iReson = 1:numel(ResonanceFields)

        % Insert new row into data arrays if the current one is not for transition iTrans!!
        %-----------------------------------------------------------------------------------
        if idxTr(iiTrans)>iTrans
          Pdat = [Pdat(1:iiTrans-1,:); ZeroRow; Pdat(iiTrans:end,:)];
          if computeIntensities
            Idat = [Idat(1:iiTrans-1,:); ZeroRow; Idat(iiTrans:end,:)];
          end
          if computeStrains
            Wdat = [Wdat(1:iiTrans-1,:); ZeroRow; Wdat(iiTrans:end,:)];
          end
          if computeGradient
            Gdat = [Gdat(1:iiTrans-1,:); ZeroRow; Gdat(iiTrans:end,:)];
          end
          idxTr = [idxTr(1:iiTrans-1); iTrans; idxTr(iiTrans:end)];
          if nPerturbNuclei>0
            for iiNuc = nPerturbNuclei:-1:1
              p_ = pPdatN{iiNuc};
              z_ = zeros(1,nOrientations,size(p_,3));
              pPdatN{iiNuc} = [p_(1:iiTrans-1,:,:); z_; p_(iiTrans:end,:,:)];
              p_ = pIdatN{iiNuc};
              pIdatN{iiNuc} = [p_(1:iiTrans-1,:,:); z_; p_(iiTrans:end,:,:)];
            end
          end
        end

        % Update position data
        %------------------------------------------------
        Pdat(iiTrans,iOri) = ResonanceFields(iReson);

        % Compute eigenvectors, eigenvalues and 1/g factor = 1/(dE/dB) if needed
        %--------------------------------------------------
        if computeEigenPairs || nPerturbNuclei>0
          % Compute resonant state vectors
          % u: lower level, v: higher level
          uv = Transitions(idxTr(iiTrans),:);

          % If eigenvectors change too much between knots, we have to
          % rediagonalize the Hamiltonian at the resonance field.
          if any(StateStability(uv,s)<Opt.RediagLimit)
            nRediags = nRediags + 1;
            if higherOrder
              [Vectors_,Energies] = gethamdata_hO(ResonanceFields(iReson),zLab_M,CoreSys,Opt.Sparse,Trans,nLevels);
              if Opt.Sparse
                [Energies,ind] = sort(diag(Energies));
                Energies = diag(Energies);
                Vectors_ = Vectors_(:,ind);
              end
            else
              if issparse(kH0)
                [Vectors_,Energies] = eigs(kH0-ResonanceFields(iReson)*kmuzL,nLevels);
                % A sort of workaround for diagonalization using eigs, the
                % energies are not ordered which results in a miscalculation
                % of mu
                [Energies,ind] = sort(diag(Energies));
                Energies = diag(Energies);
                Vectors_ = Vectors_(:,ind);
                
                %[Vectors_,Energies] = eig(full(kF+ResonanceFields(iReson)*kGzL));
              else
                [Vectors_,Energies] = eig(kH0-ResonanceFields(iReson)*kmuzL);
              end
              Energies = diag(Energies);
            end
            U = Vectors_(:,uv(1));
            V = Vectors_(:,uv(2));
            %Energies = diag(Energies);
            %V'*kGxL*U
            logmsg(3,sprintf('   %d-%d: stabilities %f and %f ===> rediagonalization',uv,StateStability(uv,s)));
          else

            z = cubicZeros(iReson);

            Ua = vectors{s}(:,uv(1)); Ub = vectors{s+1}(:,uv(1));
            [~,idx] = max(abs(Ua));
            phase = Ua(idx)/Ub(idx);
            U = Ua*(1-z) + z*phase/abs(phase)*Ub;
            U = U/norm(U);

            Va = vectors{s}(:,uv(2)); Vb = vectors{s+1}(:,uv(2));
            [~,idx] = max(abs(Va));
            phase = Va(idx)/Vb(idx);
            V = Va*(1-z) + z*phase/abs(phase)*Vb;
            V = V/norm(V);

            if computeBoltzmannPopulations
              t = cubicZeros(iReson);
              Energies = [t^3 t^2 t 1]*SplineModelCoeffs{s};
            end
          end

          if higherOrder
            if Opt.Sparse
              sp = 'sparse';
            else
              sp = '';
            end
            g1 = ham_ezho(CoreSys,[],[],sp,1);
            [g0{1},g0{2},g0{3}] = ham_ez(CoreSys,[],sp);
            if Sys.nNuclei>0
              [mu0n{1},mu0n{2},mu0n{3}] = ham_nz(CoreSys,[],sp);
              for k = 1:3
                g0{k} = g0{k} - mu0n{k};
              end
            end
            for n =3:-1:1
              kmuM{n} = -(g1{1}{n}+g0{n});
            end
            % calculate lab-frame components
            kmuzL = zLab_M(1)*kmuM{1} + zLab_M(2)*kmuM{2} + zLab_M(3)*kmuM{3};
            kmuxL = xLab_M(1)*kmuM{1} + xLab_M(2)*kmuM{2} + xLab_M(3)*kmuM{3};
            kmuyL = yLab_M(1)*kmuM{1} + yLab_M(2)*kmuM{2} + yLab_M(3)*kmuM{3};
          end
          
          % Compute dB/dE
          % dBdE is the general form of the famous 1/g factor
          % dBdE = (d(Ev-Eu)/dB)^(-1) = 1/(<v|dH/dB|v>-<u|dH/dB|u>)
          if computeFreq2Field
            dBdE = 1/abs(real((V-U)'*(-kmuzL)*(V+U)));
            % It might be quicker to take it from the first derivative
            % of the transition energy!
            %dBdE = dB(s)/abs(Diff1(iReson));
            %dBdE2 = dB(s)^2/abs(Diff2(iReson)); % second derivative
            %dBdE/dBdEold-1
            % Guard against d(Ev-Eu)/dB==0
            if dBdE>1e5
              error('1/g factor diverges because d(Ev-Eu)/dB is almost zero for transition between levels u=%d and v=%d.',uv(1),uv(2));
            end
          else
            dBdE = 1;
          end
        end

        % Calculate intensity if requested
        %--------------------------------------------------
        if computeIntensities
          
          % Compute quantum-mechanical transition rate
          mu_L = [V'*kmuxL*U; V'*kmuyL*U; V'*kmuzL*U]; % magnetic transition dipole moment, in lab frame
          if averageOverChi
            if mwmode.linearpolarizedMode
              TransitionRate = ((1-xi1^2)*norm(mu_L)^2+(3*xi1^2-1)*abs(nB0_L.'*mu_L)^2)/2;
            elseif mwmode.unpolarizedMode
              TransitionRate = ((1+xik^2)*norm(mu_L)^2+(1-3*xik^2)*abs(nB0_L.'*mu_L)^2)/4;
            elseif mwmode.circpolarizedMode
              TransitionRate = ((1+xik^2)*norm(mu_L)^2+(1-3*xik^2)*abs(nB0_L.'*mu_L)^2)/2 - ...
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
          if abs(TransitionRate)<1e-10
            TransitionRate = 0;
          end
          
          % Compute polarizations if temperature or zero-field populations are given.
          if computeBoltzmannPopulations
            Populations = exp(-BoltzmannPreFactor*(Energies-Energies(1)));
            Populations = Populations/sum(Populations);
            Polarization = Populations(u(iTrans)) - Populations(v(iTrans));
            if Polarization<0
              error('Negative thermal polarization for transition %d<->%d: %f',u(iTrans),v(iTrans),Polarization);
            end
            if nPerturbNuclei>0
              Polarization = Polarization/prod(2*Sys.I+1);            
            end
          elseif computeNonEquiPops
            switch initStateBasis
              case 'eigen'
                PopulationU = initState(uv(1),uv(1)); % lower level
                PopulationV = initState(uv(2),uv(2)); % upper level
              otherwise
                PopulationU = U'*initState*U; % lower level
                PopulationV = V'*initState*V; % upper level
            end
            Polarization = real(PopulationU - PopulationV);
            if nPerturbNuclei>0
              Polarization = Polarization/prod(2*Sys.I+1);            
            end
          else
            % no temperature given
            Polarization = 1; % same polarization for each electron transition
            Polarization = Polarization/prod(2*Sys.I+1);
          end
          
          % Update intensity results array
          Idat(iiTrans,iOri) = dBdE * TransitionRate * Polarization * photoWeight;
          % dBdE proportionality not valid near looping field coalescences!
        end
        
        % Calculate gradient of resonance frequency
        %---------------------------------------------------
        if computeGradient
          Gradient2 = real((V'-U')*kmuxL*(V+U)).^2 + real((V'-U')*kmuyL*(V+U)).^2;
          % dBdE proportionality not valid near looping field coalescences
          Gdat(iiTrans,iOri) = dBdE * ResonanceFields(iReson) * sqrt(Gradient2);
        end
        
        % Calculate width if requested
        %--------------------------------------------------
        if computeStrains
          %m = @(Op) real(V'*Op*V) - real(U'*Op*U);
          m = @(Op) real((V'-U')*Op*(V+U)); % equivalent to prev. line

          % H strain
          LineWidth2 = LineWidthSquared;
          
          % D strain
          if useDStrain
            for iEl = 1:CoreSys.nElectrons
              LineWidth2 = LineWidth2 + abs(m(dHdD{iEl}))^2;
              LineWidth2 = LineWidth2 + abs(m(dHdE{iEl}))^2;
            end
          end
          
          % g and A strain
          if usegStrain || useAStrain
            if simplegStrain
              gA2 = gAslw2{1}(:,:,iTrans);
            else
              gA2 = 0;
              for iEl = 1:Sys.nElectrons
                gA2 = gA2 + abs(m(kSzL{iEl}))*gAslw2{iEl}(:,:,iTrans);
              end
            end
            LineWidth2 = LineWidth2 + zLab_M.'*gA2*zLab_M;
          end
          
          % Convert to field value and save
          % (dBdE proportionality not valid near looping field coalescences!)
          Wdat(iiTrans,iOri) = dBdE * sqrt(LineWidth2);
        end
        %--------------------------------------------------
        
        
        % First-order approximation for nuclei
        %-------------------------------------------------------
        if nPerturbNuclei>0
          % Compute S vector expectation values for all electron spins
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
              Hc = Hquad{iiNuc} + ResonanceFields(iReson)*...
                (zLab_M(1)*Hzeem(iiNuc).x + zLab_M(2)*Hzeem(iiNuc).y + zLab_M(3)*Hzeem(iiNuc).z);
              Hu = Hu + Hc;
              Hv = Hv + Hc;
            end
            % hermitianize (important, otherwise eig returns unsorted values)
            Hu = (Hu+Hu')/2; 
            Hv = (Hv+Hv')/2;
            [Vu,dEu] = eig(Hu); dEu = diag(dEu);
            [Vv,dEv] = eig(Hv); dEv = diag(dEv);
            NucTransitionRates = abs(Vu'*Vv).^2; % the famous Mims matrix M
            % Compute and store all resonance field shifts and amplitude factors.
            % Intensity thresholds are applied later.
            [vidx,uidx] = find(ones(size(NucTransitionRates)));
            pPdatN{iiNuc}(iiTrans,iOri,:) = -dBdE*(dEv(vidx(:)) - dEu(uidx(:)));
            pIdatN{iiNuc}(iiTrans,iOri,:) = NucTransitionRates(:);
          end
        end
        %-------------------------------------------------------

        iiTrans = iiTrans + 1;
        
      end % for all resonance fields of a given transitions
      
      if ~LoopFields, break; end % only one resfield per transition ->
        % quit segment loop if resfield processed!
      
    end  % for all segments
  end % for all transitions
end % for all orientations

clear fH1 fVu fVv Hu Hv pVu pVv NucTransitionRates vidx uidx
idxTr(end) = [];
%=======================================================================

logmsg(2,'  ## %2d resonances total from %d level pairs',size(Pdat,1),nTransitions);

% Transitions post-selection
%=======================================================================

% Remove resonances with max intensity below post-selection threshold
%-----------------------------------------------------------------------
idxWeakResonances = [];
if computeIntensities
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

% Remove out-of-range resonances
%----------------------------------------------
idxOutOfRangeResonances = find(all(isnan(Pdat),2));
if ~isempty(idxOutOfRangeResonances)
  logmsg(2,'  ## %2d resonances completely out of range',numel(idxOutOfRangeResonances));
end

% Remove resonances
%----------------------------------------------
idxRmv = [idxOutOfRangeResonances; idxWeakResonances];
if numel(idxRmv)>0
  logmsg(2,'  ## removing %2d resonances (below threshold, out of range)',numel(idxRmv));
  idxTr(idxRmv) = [];
  Pdat(idxRmv,:) = [];
  if computeIntensities, Idat(idxRmv,:) = []; end
  if computeStrains, Wdat(idxRmv,:) = []; end
  if computeGradient, Gdat(idxRmv,:) = []; end
  if nPerturbNuclei>0
    for iiNuc = 1:nPerturbNuclei
      pPdatN{iiNuc}(idxRmv,:,:) = [];
    end
    if computeIntensities
      for iiNuc = 1:nPerturbNuclei
        pIdatN{iiNuc}(idxRmv,:,:) = [];
      end
    end
  end
end
Transitions = Transitions(idxTr,:);
nTransitions = size(Pdat,1);
logmsg(2,'  ## %2d resonances left',nTransitions);

if EasySpinLogLevel>=2
  partlyNaN = any(isnan(Pdat),2);
  x = full(sparse(Transitions(:,1),Transitions(:,2),double(partlyNaN)));
  nLoopPairs = sum(sum(fix(x/2)));
  nChopped = sum(partlyNaN) - nLoopPairs*2;
  if nLoopPairs>0
    logmsg(2,'  ## %2d transitions form %d looping pairs',nLoopPairs*2,nLoopPairs);
  end
  if nChopped>0
    logmsg(2,'  ## %2d transitions partly out of range',nChopped);
  end
  logmsg(2,'  ## %2d transitions fully in range',nTransitions-sum(partlyNaN));
end
logmsg(1,'  %d significant transitions with resonances in range',nTransitions);

if nTransitions==0
  fprintf(['WARNING: No resonance fields at %g GHz between %g and %g mT.\n'...
      '         Check field range and spectrometer frequency.\n'],...
    Exp.mwFreq,Exp.Range(1),Exp.Range(2));
  Output = {[],[],[],[],[]};
  varargout = Output(1:max(nargout,1));
  return
end

% Assert positive intensities, but only for thermal equilibrium populations
if ~computeNonEquiPops
  if any(Idat(:)<0)
    logmsg(-Inf,'*********** Negative intensity encountered in resfields!! Please report! **********');
  end
end
% Assert positive widths
if any(Wdat(:)<0)
  logmsg(-Inf,'*********** Negative width encountered in resfields!! Please report! **************');
end

if nMaxSegmentsReached>0
  logmsg(1,'  *** Maximum number of segments (%d) reached for %d orientations!',Opt.maxSegments,nMaxSegmentsReached);
end

% Transition post-selection for perturbational nuclei
%---------------------------------------------------------------------
if nPerturbNuclei>0
  logmsg(1,'  transition post-selection on perturbation nuclei');

  % Remove low-intensity lines from splittings
  for iiNuc = nPerturbNuclei:-1:1
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
  for iiNuc = nPerturbNuclei:-1:2
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
  if numel(Wdat)>0
    Wdat = rearrange(expand(Wdat));
  end
  if numel(Gdat)>0
    Gdat = rearrange(expand(Gdat));
  end
  Transitions = rearrange(expand(Transitions));
  clear pPdat pIdat
end

% Performance analysis
%---------------------------------------------------------------------
nResonances = nnz(~isnan(Pdat));
OldWarningState = warning;
warning off;
logmsg(2,'  ## diags %d  resonances %d  oris %d',nDiagonalizations,nResonances,nOrientations);
logmsg(1,'  ## diags/field %2.4g  diags/ori %2.4g  fields/ori %2.4g',nDiagonalizations/nResonances,...
  nDiagonalizations/nOrientations,nResonances/nOrientations);
logmsg(2,'  ## Hamiltonian rediagonalised for %d resonances (%0.3f%%)',nRediags,nRediags/nResonances*100);
logmsg(2,'  ## minimum state stability encountered %f',minStateStability);
warning(OldWarningState);

% Resonance data summary
%---------------------------------------------------------------------
logmsg(2,'  ## resonances min %g mT, max %g mT',min(Pdat(:)),max(Pdat(:)));
logmsg(2,'  ## amplitudes min %g, max %g',min(Idat(:)),max(Idat(:)));
if numel(Wdat)>0
  logmsg(2,'  ## widths min %g mT, max %g mT',min(Wdat(:)),max(Wdat(:)));
end
if numel(Gdat)>0
  logmsg(2,'  ## gradients min %g, max %g',min(Gdat(:)),max(Gdat(:)));
end

% Reshape arrays in the case of crystals with site splitting
d = dbstack;
pepperCall = numel(d)>2 && strcmp(d(2).name,'pepper');
if nSites>1 && ~pepperCall
  % Pdat, Idat, Wdat have size [nTransitions, nOrientations*nSites]
  % Resize to [nTransitions*nSites, nOrientations]
  siz = [nTransitions*nSites, numel(Pdat)/nTransitions/nSites];
  Pdat = reshape(Pdat,siz);
  if ~isempty(Idat), Idat = reshape(Idat,siz); end
  if ~isempty(Wdat), Wdat = reshape(Wdat,siz); end
  if ~isempty(Gdat), Gdat = reshape(Gdat,siz); end
end

% Sort transitions lexicograpically (unless there is more than one site)
if nSites==1
  [Transitions, I] = sortrows(Transitions);
  Pdat = Pdat(I,:);
  if ~isempty(Idat), Idat = Idat(I,:); end
  if ~isempty(Wdat), Wdat = Wdat(I,:); end
  if ~isempty(Gdat), Gdat = Gdat(I,:); end
end

% Arrange the output
Output = {Pdat,Idat,Wdat,Transitions,Gdat};
varargout = Output(1:max(nargout,1));

end
