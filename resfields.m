% resfields  Compute resonance fields for cw EPR 
%
%   Pos = resfields(Sys,Par)
%   Pos = resfields(Sys,Par,Opt)
%   [Pos,Int] = resfields(...)
%   [Pos,Int,Wid] = resfields(...)
%   [Pos,Int,Wid,Trans] = resfields(...)
%
%   Computes cw EPR line positions, intensities and widths.
%
%   Input:
%   - Sys:    spin system structure
%   - Par:    experimental parameter settings
%             mwFreq:   in GHz
%             Range:    [Bmin Bmax] in mT
%             Temperature: in K, default inf
%             Mode: 'perpendicular','parallel'
%             Orientations:  orientations of the spin system in the spectrometer
%                   2xn or 3xn array containing [phi;theta{;chi}] in
%                   radians
%   - Opt:    additonal computational options
%             Transitions, Threshold, etc
%
%   Output:
%   - Pos:    line positions
%   - Int:    line intensities, possibly including gradients
%   - Wid:    line widths
%   - Trans:  list of transitions included in the computation

function varargout = resfields(System,Params,Opt)

if (nargin==0), help(mfilename); return; end

% General
%---------------------------------------------------------------------
% Assert correct Matlab version
error(chkmlver);

% Check number of input arguments.
switch (nargin)
case 0, help(mfilename); return;
case 2, Opt = [];
case 3,
otherwise
  error('Incorrect number of inputs!');
end

if (nargout<0), error('Not enough output arguments.'); end
if (nargout>5), error('Too many output arguments.'); end

if isempty(Opt)
  Opt = struct('ununsed',NaN);
end

if ~(isstruct(System) && isstruct(Params) && isstruct(Opt))
  error('SpinSystem, Parameters and Options must be structures!');
end


% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;

% Process Spin system.
%---------------------------------------------------------------------
[System,err] = validatespinsys(System);
error(err);

DefaultSystem.lw = 0;
DefaultSystem.HStrain = [0 0 0];
DefaultSystem.gStrain = [0 0 0];
DefaultSystem.AStrain = [0 0 0];
DefaultSystem.DStrain = 0;

System = adddefaults(System,DefaultSystem);

if (System.nElectrons>1)
  if any(System.gStrain) || any(System.AStrain)
    error('Cannot use D or g/A strain in spin system with more than one electron spin.');
  end
end

if any(System.gStrain) || any(System.AStrain)
  gFull = size(System.g,1)==3*numel(System.S);
  %aFull = size(System.A,1)==3*(1+sum(System.Nucs==','));
  if gFull
    error('gStrain and AStrain are not allowed when full g matrices are given!');
  end
  if any(System.DStrain)
    error('D strain and g/A strain cannot be used at the same time.');
  end
end

if any(System.DStrain(:)) && any(System.Dpa(:))
  error('D stain cannot be used with tilted D tensors.');
end

% Process Parameters.
%---------------------------------------------------------------------
if isfield(Params,'Detection')
  error('Exp.Detection is obsolete. Use Exp.Mode instead.');
end

DefaultParams.mwFreq = NaN;
DefaultParams.Range = NaN;
DefaultParams.Orientations = NaN;
DefaultParams.CenterSweep = NaN;
DefaultParams.Temperature = inf;
DefaultParams.Mode = 'perpendicular';
DefaultParams.CrystalSymmetry = '';

Params = adddefaults(Params,DefaultParams);

if isnan(Params.Orientations)
  Params.Orientations = [0;0];
  logmsg(0,'Exp.Orientations is missing, assuming [0;0].');
end

if isnan(Params.mwFreq), error('Experiment.mwFreq is missing!'); end
mwFreq = Params.mwFreq*1e3; % GHz -> MHz

if ~isnan(Params.CenterSweep)
  if ~isnan(Params.Range)
    %logmsg(0,'Using Experiment.CenterSweep and ignoring Experiment.Range.');
  end
  Params.Range = Params.CenterSweep(1) + [-1 1]*Params.CenterSweep(2)/2;
  Params.Range = max(Params.Range,0);
end

if isfield(Params,'SearchRange'), Params.Range = Params.SearchRange; end

if isnan(Params.Range), error('Experiment.Range/Exp.CenterSweep is missing!'); end
if (diff(Params.Range)<=0) | ~isfinite(Params.Range) | ~isreal(Params.Range) | any(Params.Range<0)
  error('Params.Range is not valid!');
end

% Resonator mode
if isfield(Params,'Detection')
  error('Exp.Detection is obsolete. Use Exp.Mode instead.');
end
if strcmp(Params.Mode,'perpendicular')
  ParallelMode = 0;
elseif strcmp(Params.Mode,'parallel')
  ParallelMode = 1;
else
  error('Exp.Mode must be either ''perpendicular'' or ''parallel''.');
end
logmsg(1,'  %s mode',Params.Mode);

nPop = numel(Params.Temperature);
if (nPop>1)
  ComputeBoltzmannPopulations = 0;
  ComputeNonEquiPops = 1;
  nElectronStates = prod(2*System.S+1);
  if (nPop~=nElectronStates)
    error(sprintf('Params.Temperature must either be a scalar or a %d-vector',nElectronStates));
  end
else
  ComputeNonEquiPops = 0;
  ComputeBoltzmannPopulations = isfinite(Params.Temperature);
end


% Orientations
%-----------------------------------------------------------------------
% Contains a 2xn or 3xn array of angles (in radian units). These specify
% the relative orientation between molecular and laboratory frame
% [phi;theta;chi]. If the third angle is missing, an integration of
% the signal over the third angle is done (plane average: same B0
% direction, but B1 integrated over the plane if in perpendicular
% mode). This only affects intensity computations in perpendicular
% mode.
Orientations = Params.Orientations;
[n1,n2] = size(Orientations);
if ((n2==2)||(n2==3)) && (n1~=2) && (n1~=3)
  Orientations = Orientations.';
end
[nAngles,nOrientations] = size(Orientations);
switch nAngles
  case 2
    IntegrateOverChi = 1;
    Orientations(3,end) = 0; % Entire chi row is set to 0.
  case 3
    IntegrateOverChi = 0;
  otherwise
    error('Orientations array has %d rows instead of 2 or 3.',nAngles);
end

% Add symmetry-related sites if space group symmetry is given
if ~isempty(Params.CrystalSymmetry)
  R = sitetransforms(Params.CrystalSymmetry);
  nSites  = numel(R);
  allOrientations = zeros(nOrientations*nSites,3);
  idx = 1;
  for iOri = 1:nOrientations
    xyz0 = erot(Orientations(:,iOri)).'; % xL, yL, zL along columns
    for iSite = 1:nSites
      xyz = R{iSite}*xyz0; % active rotation
      allOrientations(idx,:) = eulang(xyz.',1);
      idx = idx + 1;
    end
  end
  Orientations = allOrientations.';
  [nAngles,nOrientations] = size(Orientations);
else
  nSites = 1;
end

% Options parsing and setting.
%---------------------------------------------------------------------
% obsolete fields
ObsoleteOptions = {''};
for iOpt = 1:numel(ObsoleteOptions)
  if isfield(Opt,ObsoleteOptions{iOpt})
    error(sprintf('Options.%s is obsolete. Please remove from code!',ObsoleteOptions{iOpt}));
  end
end

if isfield(Opt,'Perturb')
  error('Options.Perturb is obsolete. Use Opt.Method=''perturb'' or Opt.Method=''hybrid'' instead.');
end

% Fields with changed format: formerly strings/documented, now {0,1}/undocumented
changedFields = {'Intensity','Gradient'};
for iFld = 1:numel(changedFields)
  if isfield(Opt,changedFields{iFld}),
    if ischar(Opt.(changedFields{iFld}))
      error(sprintf('Options.%s is obsolete. Please remove from code!',changedFields{iFld}));
    end
  end
end

% documented fields
DefaultOptions.Transitions = [];
DefaultOptions.Threshold = 1e-4;
DefaultOptions.Hybrid = 0;
DefaultOptions.HybridNuclei = [];

% undocumented fields
DefaultOptions.nTRKnots = 3;
DefaultOptions.FuzzLevel = 1e-7;
DefaultOptions.Freq2Field = 1;
DefaultOptions.MaxKnots = 2000;
DefaultOptions.ModellingAccuracy = 2e-6;
DefaultOptions.RediagLimit = 0.95;

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

IntensitySwitch = Opt.Intensity;
GradientSwitch = Opt.Gradient;

if (Opt.Freq2Field~=1)&&(Opt.Freq2Field~=0)
  error('Options.Freq2Field incorrect!');
end
ComputeFreq2Field = Opt.Freq2Field;

ComputeStrains = (nargout>2) & ...
  any([System.HStrain(:); System.DStrain(:); System.gStrain(:); System.AStrain(:)]);
ComputeGradient = (ComputeStrains | (nargout>4)) & GradientSwitch;
ComputeIntensities = ((nargout>1) & IntensitySwitch) | ComputeGradient;
ComputeEigenPairs = ComputeIntensities | ComputeGradient | ComputeStrains;


% Preparing kernel and perturbing system Hamiltonians.
%-----------------------------------------------------------------------
logmsg(1,'- Preparations');

% :KLUDGE: Add some fuzz to the hyperfine couplings to avoid degeneracies
% if several (equivalent) nuclei are specified.
if (System.nNuclei>1)
  System.A = System.A.*(1 + Opt.FuzzLevel*rand(size(System.A)));
end

CoreSys = System;

% HFI splitting at zero field relative to mw frequency
HFIStrength = 0;
if (System.nNuclei>0)
  if System.fullA
    for iNuc = System.nNuclei:-1:1
      maxHF(iNuc) = max(max(System.A((iNuc-1)*3+(1:3))));
    end
  else
    maxHF = max(abs(System.A),[],2);
  end
  HFIStrength = maxHF(:).'.*(1/2+nucspin(System.Nucs))/mwFreq;
end

% Perturbational treatment of SHF nuclei
if (CoreSys.nNuclei>=1) && Opt.Hybrid

  %if iscell(CoreSys.Nucs), CoreSys.Nucs = CoreSys.Nucs{1}; end
  
  Nucs = nucstringparse(CoreSys.Nucs);
  
  if any(Opt.HybridNuclei>CoreSys.nNuclei)
    error('Opt.HybridNuclei is incorrect!');
  end
  perturbNuclei = ones(1,CoreSys.nNuclei);
  perturbNuclei(Opt.HybridNuclei) = 0;
  
  idx = find(perturbNuclei);
  %idx = idx & (HFIStrength<Opt.HybridHFIThreshold);
  % :TODO: Allow 1st-order PT only if (2nd-order) error smaller than field increment.  
  idx = find(idx);
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
    I = System.I(iNuc);
    nPerturbTransitions(iiNuc) = (2*I+1)^2;
    
    % Hyperfine interaction
    for iElectron = 1:System.nElectrons
      idxE = 3*(iElectron-1)+(1:3);
      if System.fullA
        A = System.A(3*(iNuc-1)+(1:3),idxE);
      else
        A = System.A(iNuc,idxE); R = eye(3);
        if isfield(System,'Apa'), R = erot(System.Apa(iNuc,idxE)); end
        A = R*diag(A)*R';
      end
      Hhfi(iElectron,iiNuc).x = A(1,1)*sop(I,1,1) + A(1,2)*sop(I,1,2) + A(1,3)*sop(I,1,3);
      Hhfi(iElectron,iiNuc).y = A(2,1)*sop(I,1,1) + A(2,2)*sop(I,1,2) + A(2,3)*sop(I,1,3);
      Hhfi(iElectron,iiNuc).z = A(3,1)*sop(I,1,1) + A(3,2)*sop(I,1,2) + A(3,3)*sop(I,1,3);
    end
    
    if ~Opt.HybridOnlyHFI
      % Nuclear Zeeman interaction
      prefactor = -nmagn/planck/1e9*System.gn(iNuc);
      Hzeem(iiNuc).x = prefactor*sop(I,1,1);
      Hzeem(iiNuc).y = prefactor*sop(I,1,2);
      Hzeem(iiNuc).z = prefactor*sop(I,1,3);
      % Nuclear quadrupole interaction
      Hquad{iiNuc} = 0;
      if (I>=1)
        Q = [0 0 0]; R = eye(3);
        if isfield(System,'Q'), Q = System.Q(iNuc,:); end
        if isfield(System,'Qpa'), R = erot(System.Qpa(iNuc,:)); end
        Q = R*diag(Q)*R';
        for c1=1:3
          for c2=1:3
            Hquad{iiNuc} = Hquad{iiNuc} + Q(c1,c2)*sop(I,1,c1)*sop(I,1,c2);
          end
        end
      end
    end

  end
  
  % Components of S vectors for computing <u|S|u>
  for iEl = System.nElectrons:-1:1
    S(iEl).x = sop(CoreSys,iEl,1);
    S(iEl).y = sop(CoreSys,iEl,2);
    S(iEl).z = sop(CoreSys,iEl,3);
  end

else
  nPerturbNuclei = 0;
end

% Hamiltonian components for the core system.
[kF,kGxM,kGyM,kGzM] = sham(CoreSys);
nCore = length(kF);
nFull = hsdim(System);
nSHFNucStates = nFull/nCore;

% Population vector for the core system
if (ComputeNonEquiPops)
  ZeroFieldPops = Params.Temperature(:);
  ZeroFieldPops = ZeroFieldPops/sum(ZeroFieldPops);
  nElStates = prod(2*System.S+1);
  ZeroFieldPops = kron(ZeroFieldPops,ones(nCore/nElStates,1));
end

if (nPerturbNuclei>0)
  logmsg(1,'  core system with %d spins and %d states',numel(spinvec(CoreSys)),nCore);
  logmsg(1,'  first-order perturbational nuclei with %d states',nSHFNucStates');
else
  if (CoreSys.nNuclei>0)
    logmsg(1,'  full treatment of all nuclei');
  end
end

% For polarized systems, pre-compute ZF eigenstates.
if (ComputeNonEquiPops)
  [ZFStates,ZFEnergies] = eig(kF);
  [ZFEnergies,idx] = sort(real(diag(ZFEnergies)));
  ZFStates = ZFStates(:,idx);
  % Correct zero-field states for S=1 and axial D
  if (CoreSys.S==1)
    if (ZFEnergies(2)==ZFEnergies(3))
      logmsg(1,'  >>>> manual zero-field states (D>0)');
      v1 = ZFStates(:,2);
      v2 = ZFStates(:,3);
      ZFStates(:,2) = (v1-v2)/sqrt(2);
      ZFStates(:,3) = (v1+v2)/sqrt(2);
    elseif (ZFEnergies(2)==ZFEnergies(1))
      logmsg(1,'  >>>> manual zero-field states (D<0)');
      v1 = ZFStates(:,1);
      v2 = ZFStates(:,2);
      ZFStates(:,2) = (v1-v2)/sqrt(2);
      ZFStates(:,1) = (v1+v2)/sqrt(2);
    end
  end
else
  ZFEnergies = sort(real(eig(kF)));
end

% Check whether looping transitions are possible.
MaxZeroFieldSplit = ZFEnergies(end) - ZFEnergies(1);
logmsg(2,'  ## maximum splitting at zero field: %g MHz',MaxZeroFieldSplit);
LoopFields = MaxZeroFieldSplit/mwFreq > 1 - 1e-6;
if (LoopFields)
  msg = '  looping transitions possible';
else
  msg = '  no looping transitions possible';
end
logmsg(1,msg);


%=======================================================================
% Transition pre-selection
%=======================================================================
% Compose a list of transitions, ie level pairs for which peak data are
% to be computed. The list is either taken from a user specified list
% in Opt.Transitions, or is composed by an automatic procedure which
% selects the most intense transitions, their number being determined
% by a threshold for the relative transitions rate. The relative transition
% rate for the most intense transition is 1.

logmsg(1,'- Transition pre-selection');  

UserTransitions = ~isempty(Opt.Transitions);
if UserTransitions

  % User-specified transitions present.
  logmsg(1,'  using %d user-specified transitions',size(Opt.Transitions,1));
  % Guarantee that lower index comes first (gives later u < v).
  if size(Opt.Transitions,2)~=2
    error('Options.Transitions must be a nx2 array of energy level indices.');
  end
  Transitions = sort(Opt.Transitions,2);
  nTransitions = size(Transitions,1);

else % Automatic pre-selection
  
  nSStates = prod(2*CoreSys.S+1);
  if (Opt.Threshold(1)==0)
    logmsg(1,'  selection threshold is zero -> using all transitions');
    TransitionRates = ones(nCore);
  elseif (CoreSys.nElectrons==1) && (CoreSys.nNuclei==0)
    logmsg(1,'  one electron spin and no nuclei -> using all transitions');
    TransitionRates = ones(nCore);
  else
    if (nOrientations>1) % if powder or multiple orientations
      % Set a coarse grid, independent of the Hamiltonian symmetry
      logmsg(1,'  selection threshold %g',Opt.Threshold(1));
      logmsg(2,'  ## (selection threshold %g, %d knots)',Opt.Threshold(1),Opt.nTRKnots);
      [phi,theta,TRWeights] = sphgrid('D2h',Opt.nTRKnots);
    else % single orientation
      phi = Params.Orientations(1);
      theta = Params.Orientations(2);
      TRWeights = 1;
    end
    % Pre-compute trigonometric functions.
    stp = sin([theta;phi].');
    ctp = cos([theta;phi].');
    centreB = mean(Params.Range); % take field at centre of scan range
    % Pre-allocate the transition rate matrix.
    TransitionRates = zeros(nCore);
    % Detector operator for transition selection.
    ExM = kGxM; EyM = kGyM; EzM = kGzM;
    % Calculate transition rates over all orientations (fixed field!).
    for iOri = 1:numel(theta)
      % Determine orientation dependent operators.
      kGpM = ctp(iOri,2)*kGxM + stp(iOri,2)*kGyM;
      EpM = ctp(iOri,2)*ExM + stp(iOri,2)*EyM;
      % Solve eigenproblem.
      [Vs,E] = eig(kF + centreB*(stp(iOri,1)*kGpM + ctp(iOri,1)*kGzM));
      % Sum up transition rates. Or take the maximum.
      if (ParallelMode)
        %TransitionRates = TransitionRates + TRWeights(iOri) * abs(Vs'*(stp(iOri,1)*EpM + ctp(iOri,1)*EzM)*Vs).^2;
        TransitionRates = max(TransitionRates,TRWeights(iOri) * abs(Vs'*(stp(iOri,1)*EpM + ctp(iOri,1)*EzM)*Vs).^2);
      else % perpendicular
        EyL = -stp(iOri,2)*ExM + ctp(iOri,2)*EyM;
        ExL =  ctp(iOri,1)*EpM - stp(iOri,1)*EzM;
        %TransitionRates = TransitionRates + TRWeights(iOri) * (abs(Vs'*ExL*Vs).^2 + abs(Vs'*EyL*Vs).^2);
        TransitionRates = max(TransitionRates,TRWeights(iOri) * (abs(Vs'*ExL*Vs).^2 + abs(Vs'*EyL*Vs).^2));
      end
      %Vectors{iOri} = Vs;
    end
    % Free unused memory.
    clear Vs E idx kGpM ExM EyM EzM EpM; % Vectors
  end
  
  % Remove lower triangular part
  idxLowerTriangle = logical(tril(ones(nCore)));
  % Remove nuclear transitions
  if (max(HFIStrength)<0.5) && (Opt.Threshold(1)>0)
    idxNuclearTransitions = logical(kron(eye(nSStates),ones(nCore/nSStates)));
    keepidx = ~(idxLowerTriangle | idxNuclearTransitions);
  else
    keepidx = ~(idxLowerTriangle);
  end
  TransitionRates(~keepidx) = [];
  [u,v] = find(keepidx); % Compute level pair indices.
  Transitions = [u,v];
  clear keepidx u v idxLowerTriangle idxNuclearTransitions;

  % Use threshold for number determination.
  nTransitions = sum(TransitionRates>Opt.Threshold(1)*max(TransitionRates));
  
  % Sort TransitionRates in descending order!
  [unused,idx] = sort(-TransitionRates);
  % Select most intense transitions.
  Transitions = Transitions(idx(1:nTransitions),:);
  clear unused TransitionRates idx;
  
end

% Terminate if the transition list is empty.
if isempty(Transitions)
  error('No transitions selected! Decrease Opt.Threshold.');
end

if any(Transitions(:)>nCore)
  error('Level index in Options.Transitions is out of range.');
end

% Diagnostic display.
logmsg(1,'  %d transitions pre-selected',nTransitions);

% Compute indices and variables used later in the algorithm
u = Transitions(:,1);
v = Transitions(:,2); % u < v
nTransitions = length(u);
upTRidx = u + (v-1)*nCore; % Indices into UPPER triangle.
%loTRidx = v + (u-1)*nCore; % Indices into LOWER triangle.
Trans = upTRidx; % One-number transition indices.

% Now, if the transitions were selected automatically, they are in
% descending order according to their average intensity. If user-
% specified, their order has not been changed.
%=======================================================================


%=======================================================================
% Line width preparations
%=======================================================================
logmsg(1,'- Broadenings',nTransitions);
if (ComputeStrains)
  logmsg(1,'  using strains',nTransitions);
  
  % Frequency domain residual width tensor.
  %-----------------------------------------------
  HStrain2 = CoreSys.HStrain.^2;
  
  % D strain
  %-----------------------------------------------
  % Diagonal of D tensor: [-D/3+E, -D/3-E, 2D/3] in MHz.
  % H = D*(Sz^2-S(S+1)/3) + E*(Sx^2-Sy^2)
  %   = D*(2*Sz^2-Sx^2-Sy^2)/3 + E*(Sx^2-Sy^2)
  % D strain: independent distributions in D and E parameters.
  % System.DStrain: [FWHM_D FWHM_E] in MHz.
  UseDStrain = any(CoreSys.DStrain(:));
  if UseDStrain
    for iEl = 1:CoreSys.nElectrons
      
      % Since S^T D S = S^T (R D_diag R^T) S = (S^T R) D_diag (R^T S),
      % we have to compute R^T S = (SxD SyD SzD)^T to stay in the
      % eigenframe of D for the D and E strain.
      R = eye(3);
      if any(CoreSys.Dpa(iEl,:))
        R = erot(CoreSys.Dpa(iEl,:));
      end
      R = R';
      
      % Construct Zeeman basis operators
      Sx_ = sop(CoreSys,iEl,1);
      Sy_ = sop(CoreSys,iEl,2);
      Sz_ = sop(CoreSys,iEl,3);

      % Compute squared cartesian spin operators in D eigenframe
      SxD2_ = (R(1,1)*Sx_ + R(1,2)*Sy_ + R(1,3)*Sz_)^2;
      SyD2_ = (R(2,1)*Sx_ + R(2,2)*Sy_ + R(2,3)*Sz_)^2;
      SzD2_ = (R(3,1)*Sx_ + R(3,2)*Sy_ + R(3,3)*Sz_)^2;
      
      % Tranform correlated D-E strain to uncorrelated coordinates F and G
      r = 0; % correlation coefficient between D and E
      if size(CoreSys.DStrain,2)==3
        r = CoreSys.DStrain(iEl,3);
      end
      DeltaD = CoreSys.DStrain(iEl,1);
      DeltaE = CoreSys.DStrain(iEl,2);
      if (r~=0)
        logmsg(1,'  correlated D strain for electron spin %d',iEl);
        R11 = DeltaD^2;
        R22 = DeltaE^2;
        R12 = r*DeltaD*DeltaE;
        R = [R11 R12; R12 R22]; % covariance matrix
        [V,L] = eig(R); L = sqrt(diag(L));
        DeltaF = L(1);
        DeltaG = L(2);
        F_  = V(1,1)*(2*SzD2_-SxD2_-SyD2_)/3  + V(1,2)*(SxD2_-SyD2_);
        G_  = V(2,1)*(2*SzD2_-SxD2_-SyD2_)/3  + V(2,2)*(SxD2_-SyD2_);
      else
        DeltaF = DeltaD;
        DeltaG = DeltaE;
        F_  = (2*SzD2_-SxD2_-SyD2_)/3;
        G_  = (SxD2_-SyD2_);
      end
      % Compute Hamiltonian derivatives, pre-multiply with
      % strain FWHMs.
      dHdD{iEl} = DeltaF * F_;
      dHdE{iEl} = DeltaG * G_;

    end
    clear Sx_ Sy_ Sz_ SxD2_ SyD2_ SzD2_ F_ G_;
  end
  
  % g-A strain
  %-------------------------------------------------
  % g strain tensor is taken to be along the g tensor itself.
  gslw = diag(CoreSys.gStrain./CoreSys.g(1,:)*mwFreq);
  if (CoreSys.nNuclei>0) && any(CoreSys.AStrain)
    % Get A strain matrix in g frame.
    AStrain = diag(CoreSys.AStrain);
    if isfield(CoreSys,'Apa')
      Rp = erot(CoreSys.Apa(1,:));
      AStrain = Rp*AStrain*Rp.';
    end
    % Diagonalize Hamiltonian at centre field.
    centreB = mean(Params.Range);
    [Vs,E] = eig(kF + centreB*kGzM);
    [E,idx] = sort(real(diag(E)));
    Vs = Vs(:,idx);
    % Calculate effective mI of nucleus 1 for all eigenstates.
    mI = real(diag(Vs'*sop(CoreSys,2,3)*Vs));
    mITr = mean(mI(Transitions),2);
    % compute A strain array
    Aslw = reshape(mITr(:,ones(1,9)).',[3,3,nTransitions]).*...
      repmat(AStrain,[1,1,nTransitions]);
    gAslw2 = (repmat(gslw,[1,1,nTransitions])+Aslw).^2;
    clear Aslw Vs E idx mI mITr
  else
    gAslw2 = repmat(gslw.^2,[1,1,nTransitions]);
  end
  clear gslw
  % gAslw2 = now an 3D array with 3x3 strain line-width matrices
  % for each transition piled up along the third dimension.
  gA = any(gAslw2(:)); % Switch to indicate g/A strain.
  
  if any(HStrain2), logmsg(2,'  ## using H strain'); end
  if gA, logmsg(2,'  ## using g/A strain'); end
  if UseDStrain, logmsg(2,'  ## using D strain'); end
  
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
Pdat = ones(nTransitions,nOrientations)*NaN;
if (nPerturbNuclei>0)
  for iiNuc = nPerturbNuclei:-1:1
    pPdatN{iiNuc} = zeros(nTransitions,nOrientations,nPerturbTransitions(iiNuc));
  end
end

if (ComputeIntensities)
  Idat = ones(nTransitions,nOrientations)*NaN;
  if (nPerturbNuclei>0)
    %qqq pIdatN = cell(nTransitions,nOrientations,nPerturbNuclei);
    for iiNuc = nPerturbNuclei:-1:1
      pIdatN{iiNuc} = zeros(nTransitions,nOrientations,nPerturbTransitions(iiNuc));
    end
  end
  if ComputeBoltzmannPopulations
    % Pre-factor for thermal equilibrium populations computations.
    BoltzmannPreFactor = -1e6*planck/boltzm/Params.Temperature; % MHz^-1
  end
  msg = [msg ', intensities'];
else
  Idat = ones(nTransitions,nOrientations);
end

if (ComputeStrains)
  msg = [msg ', strain widths'];
  Wdat = ones(nTransitions,nOrientations)*NaN;
else
  Wdat = [];
end

% Pre-allocate Gradient if it is requested.
if (ComputeGradient)
  msg = [msg ', gradients'];
  Gdat = zeros(nTransitions,nOrientations);
else
  Gdat = [];
end
logmsg(1,'- Resonance data: computing %s',msg);

idxTr = [(1:nTransitions).'; inf];

% Preparation for the adaptive iterative bisection
%---------------------------------------------------------
Accuracy = mwFreq*Opt.ModellingAccuracy;

M = [2 -2 1 1; -3 3 -2 -1; 0 0 1 0; 1 0 0 0];
ZeroRow = NaN*ones(1,nOrientations);
maxSlope = max(max([eig(kGxM) eig(kGyM) eig(kGzM)]));
nDiagonalizations = 0; % or 4? (1 for F and 3 for maxSlope)
nRediags = 0;
minStateStability = 2;
nMaxKnotsReached = 0;

logmsg(2,'  ## modelling accuracy %0.3f MHz, max slope %g MHz/mT',Accuracy,maxSlope);

% Loop over all given orientations.
%-----------------------------------------------------------------------

startTime = cputime;
logstr = '';
for iOri = 1:nOrientations

  if (EasySpinLogLevel>=1)
    if (iOri>1)
      remainingTime = (cputime-startTime)/(iOri-1)*(nOrientations-iOri+1);
      backspace = repmat(sprintf('\b'),1,numel(logstr));
      hours = fix(remainingTime/3600);
      minutes = fix(remainingTime/60 - 60*hours);
      seconds = remainingTime - 3600*hours - 60*minutes;
      logstr = sprintf('  %d/%d orientations, remaining time %02d:%02d:%0.1f\n', ...
        iOri, nOrientations, hours, minutes, seconds);
      fprintf([backspace logstr]);
    else
      if (nOrientations>1)
        logstr = sprintf('  1/%d orientations, remaining time unknown\n',nOrientations);
        fprintf(logstr);
      end
    end
  end
  
  % Set up Hamiltonians for 3 lab principle orientations
  %-----------------------------------------------------
  [xLab,yLab,zLab] = erot(Orientations(:,iOri));
  
  % z laboratoy axis: external static field
  kGzL = zLab(1)*kGxM + zLab(2)*kGyM + zLab(3)*kGzM;
  % x laboratory axis: mw excitation field
  kGxL = xLab(1)*kGxM + xLab(2)*kGyM + xLab(3)*kGzM;
  % y laboratory vector: needed for gradient calculation
  % and the integration over all mw field orientations.
  kGyL = yLab(1)*kGxM + yLab(2)*kGyM + yLab(3)*kGzM;

  if ComputeStrains
    LineWidthSquared = zLab.^2*HStrain2.';
  end
  
  %===========================================================
  % Iterative bisection
  %-----------------------------------------------------------  

  % V: eigenvectors, E: energies, D1: dE/dB
  % dE: transition energies for transitions in list Trans
  clear VV E D1 dE;
  Knots = Params.Range;
  [VV{2},E{2},D1{2},dE{2}] = gethamdata(Knots(2),kF,kGzL,Trans);
  [VV{1},E{1},D1{1},dE{1}] = gethamdata(Knots(1),kF,kGzL,Trans);
  nDiagonalizations = nDiagonalizations+2;
  UnFinished = 1;

  % Loop until maximum allowed number of knots is reached.
  while any(UnFinished) && (numel(Knots)<Opt.MaxKnots)
    
    s = find(UnFinished); s = s(1); % Find unfinished segment
    dB = Knots(s+1)-Knots(s);
  
    if (E{s+1}(end)-E{s+1}(1) > mwFreq)
      if (LoopFields)
        ResonancePossible = abs((dE{s}+dE{s+1})/2-mwFreq) <= maxSlope*dB;
      else
        ResonancePossible = (dE{s}-mwFreq).*(dE{s+1}-mwFreq) <= 0;
      end
    else
      ResonancePossible = 0;
    end
    
    if any(ResonancePossible)
      % compute error
      newB = (Knots(s)+Knots(s+1))/2;
      [Ve,En,Di1,dEn] = gethamdata(newB,kF,kGzL,Trans);
      nDiagonalizations = nDiagonalizations+1;
      Error = 2*(1/2*(E{s}+E{s+1}) + dB/8*(D1{s}-D1{s+1}) - En);
      Incl = zeros(1,nCore);
      Incl([u(ResonancePossible) v(ResonancePossible)]) = 1;
      Error = abs(Error(logical(Incl)));
      % bisect
      Knots = [Knots(1:s) newB Knots(s+1:end)];
      E = [E(1:s) {En} E(s+1:end)];
      VV = [VV(1:s) {Ve} VV(s+1:end)];
      D1 = [D1(1:s) {Di1} D1(s+1:end)];
      dE = [dE(1:s) {dEn} dE(s+1:end)];
      % mark unfinished or finished depending on error
      UnFinished = [UnFinished(1:s-1) (max(Error)>Accuracy)*[1 1] UnFinished(s+1:end)];
    else
      UnFinished(s) = 0; % mark segment as finished
    end
    
  end
  
  %------------------------------------------------------------
  
  if numel(Knots)>=Opt.MaxKnots
    nMaxKnotsReached = nMaxKnotsReached + 1;
  end
  
  nKnots = numel(Knots);
  nSegs = nKnots-1;
  dB = diff(Knots);

  logmsg(2,'   segmentation finished, %d knots, %d segments',nKnots,nSegs);
  
  % Compute eigenvector cross products to determine how strongly eigenvectors
  % change over a segment.
  for s = nSegs:-1:1
    %StateStability(:,s) = abs(diag(VV{s+1}'*VV{s}));
    StateStability(:,s) = abs(sum(conj(VV{s+1}).*VV{s})).';
  end

  % Cubic polynomial coefficients of the entire spline model of the
  % energy level diagram.
  Coeffs = cell(1,nSegs);
  for s = 1:nSegs
    Coeffs{s} = M*[E{s}; E{s+1}; dB(s)*D1{s}; dB(s)*D1{s+1}];
  end
  
  %fprintf('Orientation %d/%d: %d segments\n',iOri,nOrientations,nSegs);
  iiTrans = 1;
  for iTrans = 1:nTransitions % run over all level pairs

    % Find first position of transition iTrans in idxTr.
    while idxTr(iiTrans)<iTrans, iiTrans=iiTrans+1; end

    % Loop over all field segments and compute resonance fields
    %-----------------------------------------------------------
    for s = 1:nSegs

      % Construct cubic polynomial coeffs of resonance function Ev-Eu-freq
      Co = Coeffs{s};
      C = Co(:,v(iTrans)) - Co(:,u(iTrans));
      C(4) = C(4) - mwFreq;

      % Find zeros, first and second derivatives of E(v)-E(u)-mwFreq
      %[Zeros,Diff1,Diff2] = cubicsolve(C,LoopFields);

      % Problem here: cubicsolve should only be called if a resonance
      % is possible. This can be determined as above using maxSlope etc.
      % cubicsolve takes too long to find this out (esp. for LoopFields=1).
      Zeros = cubicsolve(C,LoopFields);
      if isempty(Zeros), continue, end

      ResonanceFields = Knots(s) + dB(s)*Zeros;

      for iReson=1:numel(ResonanceFields) % loop over all resonances for transition iTrans in the current segment

        % Insert new row into data arrays if the current one is not for transition iTrans!!
        %-----------------------------------------------------------------------------------
        if idxTr(iiTrans)>iTrans
          Pdat = [Pdat(1:iiTrans-1,:); ZeroRow; Pdat(iiTrans:end,:)];
          if (ComputeIntensities)
            Idat = [Idat(1:iiTrans-1,:); ZeroRow; Idat(iiTrans:end,:)];
          end
          if (ComputeStrains)
            Wdat = [Wdat(1:iiTrans-1,:); ZeroRow; Wdat(iiTrans:end,:)];
          end
          if (ComputeGradient)
            Gdat = [Gdat(1:iiTrans-1,:); ZeroRow; Gdat(iiTrans:end,:)];
          end
          idxTr = [idxTr(1:iiTrans-1); iTrans; idxTr(iiTrans:end)];
          if (nPerturbNuclei>0)
            for iiNuc = nPerturbNuclei:-1:1
              p_ = pPdatN{iiNuc};
              z_ = zeros(1,nOrientations,size(p_,3));
              pPdatN{iiNuc} = [p_(1:iiTrans-1,:,:); z_; ...
                               p_(iiTrans:end,:,:)];
              p_ = pIdatN{iiNuc};
              pIdatN{iiNuc} = [p_(1:iiTrans-1,:,:); z_; ...
                               p_(iiTrans:end,:,:)];
            end
          end
        end

        % Update position data
        %------------------------------------------------
        Pdat(iiTrans,iOri) = ResonanceFields(iReson);

        % Compute eigenvectors, eigenvalues and 1/(dE/dB) if needed.
        %--------------------------------------------------
        if ComputeEigenPairs || (nPerturbNuclei>0)
          % Compute resonant state vectors
          % u: lower level, v: higher level
          uv = Transitions(idxTr(iiTrans),:);

          % If eigenvectors change too much between knots, we have to
          % rediagonalize the Hamiltonian at the resonance field.
          if any(StateStability(uv,s)<Opt.RediagLimit)
            nRediags = nRediags + 1;
            [Vectors,Energies] = eig(kF+ResonanceFields(iReson)*kGzL);
            U = Vectors(:,uv(1));
            V = Vectors(:,uv(2));
            Energies = diag(Energies);
            %V'*kGxL*U
            logmsg(3,sprintf('   %d-%d: stabilities %f and %f ===> rediagonalization',uv,StateStability(uv,s)));
          else

            z = Zeros(iReson);

            Ua = VV{s}(:,uv(1)); Ub = VV{s+1}(:,uv(1));
            [val,idx] = max(abs(Ua));
            phase = Ua(idx)/Ub(idx);
            U = Ua*(1-z) + z*phase/abs(phase)*Ub;
            U = U/norm(U);

            Va = VV{s}(:,uv(2)); Vb = VV{s+1}(:,uv(2));
            [val,idx] = max(abs(Va));
            phase = Va(idx)/Vb(idx);
            V = Va*(1-z) + z*phase/abs(phase)*Vb;
            V = V/norm(V);

            t = Zeros(iReson);
            if ComputeBoltzmannPopulations
              t = [t^3 t^2 t 1];
              Energies = t*Coeffs{s};
            end
          end

          % Compute dB/dE
          % dBdE is the famous 1/g factor, in general
          % (d(Ev-Eu)/dB)^(-1) = 1/(<v|dH/dB|v>-<u|dH/dB\u>)
          if (ComputeFreq2Field)
            dBdE = 1/abs(real((V-U)'*kGzL*(V+U)));
            % It might be quicker to take it from the first derivative
            % of the transition energy!
            %dBdE = dB(s)/abs(Diff1(iReson));
            %dBdE2 = dB(s)^2/abs(Diff2(iReson)); % second derivative
            %dBdE/dBdEold-1
          else
            dBdE = 1;
          end
        end

        % Calculate intensity if requested.
        %--------------------------------------------------
        if (ComputeIntensities)
          
          % Compute quantum-mechanical transition rate
          if (ParallelMode)
            if (IntegrateOverChi)
              % Transition propability integrated over chi, the third
              % Euler angle between lab and mol frames. 
              TransitionRate = (2*pi)*abs(V'*kGzL*U).^2;
            else
              TransitionRate = abs(V'*kGzL*U).^2;
            end
          else % perpendicular mode
            if (IntegrateOverChi)
              % Transition propability integrated over chi, the third
              % Euler angle between lab and mol frames. 
              TransitionRate = pi * (abs(V'*kGxL*U).^2 + abs(V'*kGyL*U).^2);
            else
              % Transition propability along x lab axis, this is the line
              % intensity for a single (phi,theta,chi) orientation.
              TransitionRate = abs(V'*kGxL*U).^2;
            end
          end
          
          % Compute polarization if temperature or zero-field populations are given.
          if (ComputeBoltzmannPopulations)
            Populations = exp(BoltzmannPreFactor*(Energies-Energies(1)));
            Polarization = (Populations(u(iTrans)) - Populations(v(iTrans)))/sum(Populations);
            if (nPerturbNuclei>0)
              Polarization = Polarization/prod(2*System.I+1);            
            end
          elseif (ComputeNonEquiPops)
            PopulationU = (abs(ZFStates'*U).^2).'*ZeroFieldPops; % lower level
            PopulationV = (abs(ZFStates'*V).^2).'*ZeroFieldPops; % upper level
            Polarization = PopulationU - PopulationV;
          else
            % no temperature given
            Polarization = 1; % default polarization for each electron transition
            Polarization = Polarization/prod(2*System.I+1);
          end
          
          % Update intensity results array
          Idat(iiTrans,iOri) = dBdE * TransitionRate * Polarization;
          % dBdE proportionality not valid near looping field coalescences
        end
        
        % Calculate gradient of resonance frequency.
        %---------------------------------------------------
        if (ComputeGradient)
          Gradient2 = real((V'-U')*kGxL*(V+U)).^2 + real((V'-U')*kGyL*(V+U)).^2;
          % dBdE proportionality not valid near looping field coalescences
          Gdat(iiTrans,iOri) = dBdE * ResonanceFields(iReson) * sqrt(Gradient2);
        end
        
        % Calculate width if requested.
        %--------------------------------------------------
        if ComputeStrains
          % H strain
          LineWidth2 = LineWidthSquared;
          % g and A strain
          if gA
            LineWidth2 = LineWidth2 + zLab*gAslw2(:,:,iTrans)*zLab.';
          end
          % D strain
          if UseDStrain
            for iEl = 1:CoreSys.nElectrons
              % add D contribution
              LineWidth2 = LineWidth2 + abs(real((V'-U')*dHdD{iEl}*(V+U)))^2;
              % add E contribution
              LineWidth2 = LineWidth2 + abs(real((V'-U')*dHdE{iEl}*(V+U)))^2;
            end
          end
          % Convert to field value and save
          % (dBdE proportionality not valid near looping field
          % coalescences!)
          Wdat(iiTrans,iOri) = dBdE * sqrt(LineWidth2);
        end
        %--------------------------------------------------
        
        
        % First-order approximation for nuclei
        %-------------------------------------------------------
        if (nPerturbNuclei>0)
          % Compute S vector expectation values for all electron spins
          for iEl = System.nElectrons:-1:1
            Su(:,iEl) = [U'*S(iEl).x*U; U'*S(iEl).y*U; U'*S(iEl).z*U];
            Sv(:,iEl) = [V'*S(iEl).x*V; V'*S(iEl).y*V; V'*S(iEl).z*V];
          end
          % Build and diagonalize nuclear sub-Hamiltonians
          for iiNuc = 1:nPerturbNuclei
            Hu = 0;
            Hv = 0;
            % Hyperfine (dependent on S)
            for iEl = 1:System.nElectrons
              Hu = Hu + Su(1,iEl)*Hhfi(iEl,iiNuc).x + Su(2,iEl)*Hhfi(iEl,iiNuc).y + Su(3,iEl)*Hhfi(iEl,iiNuc).z;
              Hv = Hv + Sv(1,iEl)*Hhfi(iEl,iiNuc).x + Sv(2,iEl)*Hhfi(iEl,iiNuc).y + Sv(3,iEl)*Hhfi(iEl,iiNuc).z;
            end
            % Nuclear Zeeman and quadrupole (independent of S)
            if ~Opt.HybridOnlyHFI
              Hc = Hquad{iiNuc} + ResonanceFields(iReson)*...
                (zLab(1)*Hzeem(iiNuc).x + zLab(2)*Hzeem(iiNuc).y + zLab(3)*Hzeem(iiNuc).z);
              Hu = Hu + Hc;
              Hv = Hv + Hc;
            end
            % symmetrize (important, otherwise eig returns unsorted values)
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
      
      if (~LoopFields), break; end % only one resfield per transition ->
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
if ComputeIntensities
  if numel(Opt.Threshold)==1
    PostSelectionThreshold = Opt.Threshold(1);
  else
    PostSelectionThreshold = Opt.Threshold(2);
  end
  maxIntensity = max(abs(Idat),[],2); % maximum over orientations
  absoluteMaxIntensity = max(maxIntensity);
  idxWeakResonances = find(maxIntensity<absoluteMaxIntensity*PostSelectionThreshold);
  if ~isempty(idxWeakResonances)
    logmsg(2,'  ## %2d resonances below relative intensity threshold %f',numel(idxWeakResonances),PostSelectionThreshold);
  else
    logmsg(2,'  ## all resonances above relative intensity threshold %f',numel(idxWeakResonances),PostSelectionThreshold);
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
  Idat(idxRmv,:) = [];
  if (ComputeStrains), Wdat(idxRmv,:) = []; end
  if (ComputeGradient), Gdat(idxRmv,:) = []; end
  if (nPerturbNuclei>0)
    for iiNuc = 1:nPerturbNuclei
      pPdatN{iiNuc}(idxRmv,:,:) = [];
      pIdatN{iiNuc}(idxRmv,:,:) = [];
    end
  end
end
Transitions = Transitions(idxTr,:);
nTransitions = size(Pdat,1);
logmsg(2,'  ## %2d resonances left',nTransitions);

if (EasySpinLogLevel>=2),
  partlyNaN = any(isnan(Pdat),2);
  x = full(sparse(Transitions(:,1),Transitions(:,2),double(partlyNaN)));
  nLoopPairs = sum(sum(fix(x/2)));
  nChopped = sum(partlyNaN) - nLoopPairs*2;
  if (nLoopPairs>0)
    logmsg(2,'  ## %2d transitions form %d looping pairs',nLoopPairs*2,nLoopPairs);
  end
  if (nChopped>0)
    logmsg(2,'  ## %2d transitions partly out of range',nChopped);
  end
  logmsg(2,'  ## %2d transitions fully in range',nTransitions-sum(partlyNaN));
end
logmsg(1,'  %d significant transitions with resonances in range',nTransitions);

if (nTransitions==0)
  fprintf(['WARNING: No resonance fields at %g GHz between %g and %g mT.\n'...
      '         Check field range and spectrometer frequency.\n'],...
    Params.mwFreq,Params.Range(1),Params.Range(2));
  Output = {[],[],[],[],[]};
  varargout = Output(1:max(nargout,1));
  return
end

% Assert positive widths!
if any(Wdat(:)<0),
  logmsg(-inf,'*********** Negative width encountered in resfields!! Please report! **************');
end

% Assert positive intensities, but only for thermal equilibrium populations
if (~ComputeNonEquiPops)
  if any(Idat(:)<0),
    logmsg(-inf,'*********** Negative intensity encountered in resfields!! Please report! **********');
  end
end

if (nMaxKnotsReached>0)
  logmsg(1,'  *** Maximum number of knots (%d) reached for %d orientations!',Opt.MaxKnots,nMaxKnotsReached);
end

% Transition post-selection for perturbational nuclei
%---------------------------------------------------------------------
if (nPerturbNuclei>0)
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
if (nPerturbNuclei>0)
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
%---------------------------------------------------------------------
if (nPerturbNuclei>0)
  nSubTransitions = size(pPdat,3); % transitions per electronic level pair
  nTotalTrans = nTransitions*nSubTransitions;
  Pdat = reshape(permute(repmat(Pdat,[1,1,nSubTransitions]) +pPdat,[1 3 2]),nTotalTrans,nOrientations);
  Idat = reshape(permute(repmat(Idat,[1,1,nSubTransitions]).*pIdat,[1 3 2]),nTotalTrans,nOrientations);
  if numel(Gdat)>0
    Gdat = reshape(permute(repmat(Gdat,[1,1,nSubTransitions]),[1 3 2]),nTotalTrans,nOrientations);
  end
  if numel(Wdat)>0
    Wdat = reshape(permute(repmat(Wdat,[1,1,nSubTransitions]),[1 3 2]),nTotalTrans,nOrientations);
  end
  Transitions = reshape(permute(repmat(Transitions,[1 1 nSubTransitions]),[1 3 2]),nTotalTrans,2);
  clear pPdat pIdat
end

% Performance analysis
%---------------------------------------------------------------------
nResonances = sum(sum(~isnan(Pdat)));
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
if (nSites>1) && ~isfield(Opt,'peppercall')
  siz = [nTransitions*nSites, numel(Pdat)/nTransitions/nSites];
  Pdat = reshape(Pdat,siz);
  if ~isempty(Idat), Idat = reshape(Idat,siz); end
  if ~isempty(Wdat), Wdat = reshape(Wdat,siz); end
  if ~isempty(Gdat), Gdat = reshape(Gdat,siz); end
end

% Arrange the output.
Output = {Pdat,Idat,Wdat,Transitions,Gdat};
varargout = Output(1:max(nargout,1));

return
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


function [V,E,D1,dE] = gethamdata(B,F,G,idx)

% Compute eigenpairs of Hamiltonian
[V,E] = eig(F+B*G);
E = diag(E).';

% Compute correct eigenvectors for zero-field degeneracies
if (B==0)
  dE = abs(diff(E)).';
  tol = 1e3*eps*max(dE);
  blk = cumsum([1; dE>tol]);
  GG = V'*G*V;
  GG = (GG+GG')/2; % important: symmetrise
  VV_ = [];
  for k=1:max(blk)
    ix = find(blk==k);
    [v,e] = eig(GG(ix,ix));
    VV_ = blkdiag(VV_,v);
  end
  V = V*VV_;
end

% Compute first derivative dE/dB
D1 = real(diag(V'*G*V)).';

% Compute transition frequencies Ev-Eu
M = E(:);
M = M(:,ones(1,length(E)));
dE = M.' - M;
dE = dE(idx);

return
