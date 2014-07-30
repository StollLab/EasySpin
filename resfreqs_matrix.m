% resfreqs_matrix Compute resonance frequencies for frequency-domain EPR
%
%   Pos = resfreqs_matrix(Sys)
%   Pos = resfreqs_matrix(Sys,Exp)
%   Pos = resfreqs_matrix(Sys,Exp,Opt)
%   [Pos,Int] = resfreqs_matrix(...)
%   [Pos,Int,Wid] = resfreqs_matrix(...)
%   [Pos,Int,Wid,Trans] = resfreqs_matrix(...)
%
%   Computes frequency-domain EPR line positions, intensities and widths.
%
%   Input:
%   - Sys:    spin system structure
%   - Exp:    experimental parameter settings
%             Range:    [nu_min nu_max] in GHz
%             Field:    B_ 0 in mT, 0 by default
%             Temperature: in K, by default off (NaN)
%             Mode:         resonator mode: 'perpendicular' (default), 'parallel', [k_tilt alpah_pol]
%             Polarization: 'linear' (default), 'circular+', 'circular-', 'unpolarized'
%   - Opt:    additonal computational options
%             Transitions, Threshold, etc
%
%   Output:
%   - Pos:    line positions
%   - Int:    line intensities
%   - Wid:    line widths
%   - Trans:  list of transitions included in the computation

function varargout = resfreqs_matrix(System,Exp,Opt)

if (nargin==0), help(mfilename); return; end

% General
%---------------------------------------------------------------------
% Assert correct Matlab version
error(chkmlver);

% Guard against wrong number of input or output arguments.
if (nargin<1), error('Please supply a spin system as first parameter.'); end
if (nargin>3), error('Too many input arguments, the maximum is three.'); end
%Initialize options structure to zero if not given.
if (nargin<2), Exp = struct('unused',NaN); end
if (nargin<3), Opt = struct('unused',NaN); end
if isempty(Opt), Opt = struct('unused',NaN); end

if (nargout<0), error('Not enough output arguments.'); end
if (nargout>4), error('Too many output arguments.'); end

if isempty(Exp)
  Exp = struct('unused', NaN);
end
if isempty(Opt)
  Opt = struct('ununsed',NaN);
end

if ~(isstruct(System) && isstruct(Exp) && isstruct(Opt))
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
DefaultSystem.A=0;
DefaultSystem.HStrain = [0 0 0];
DefaultSystem.gStrain = [0 0 0];
DefaultSystem.AStrain = [0 0 0];
DefaultSystem.DStrain = 0;
DefaultSystem.gAStrainCorr = +1;

System = adddefaults(System,DefaultSystem);

if (numel(System.gAStrainCorr)~=1) || ~isnumeric(System.gAStrainCorr) || ...
    (System.gAStrainCorr==0) || ~isfinite(System.gAStrainCorr)
  error('Sys.gAStrainCorr must be a single number, either +1 or -1.');
end
System.gAStrainCorr = sign(System.gAStrainCorr);

if (System.nElectrons>1)
  if any(System.gStrain) || any(System.AStrain)
    error('Cannot use D or g/A strain in spin system with more than one electron spin.');
  end
end

if any(System.gStrain) || any(System.AStrain)
  gFull = size(System.g,1)==3*numel(System.S);
  if gFull
    error('gStrain and AStrain are not allowed when full g matrices are given!');
  end
  if any(System.DStrain)
    error('D strain and g/A strain cannot be used at the same time.');
  end
end

if any(System.DStrain(:)) && any(System.DFrame(:))
  error('D stain cannot be used with tilted D tensors.');
end


% Process Parameters.
%---------------------------------------------------------------------
DefaultExp.Range = NaN;
DefaultExp.Field = NaN;
DefaultExp.CenterSweep = NaN;
DefaultExp.Temperature = NaN;
DefaultExp.Mode = 'perpendicular';
DefaultExp.Polarization = '';

DefaultExp.CrystalOrientation = [];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];

Exp = adddefaults(Exp,DefaultExp);

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
if (diff(Exp.Range)<=0) | ~isfinite(Exp.Range) | ~isreal(Exp.Range) | any(Exp.Range<0)
  error('Exp.Range is not valid!');
end

Exp.Range = Exp.Range*1e3; % GHz -> MHz, for comparison with Pdat

% Determine excitation mode
p_excitationgeometry;

% Temperature
if ~numel(Exp.Temperature)
  error('Specify a single number in Exp.Temperature');
else
  if isinf(Exp.Temperature)
    error('If given, Exp.Temperature must be a finite value.');
  end
  ComputeBoltzmannPopulations = ~isnan(Exp.Temperature);
end


% Process crystal orientations, crystal symmetry, and frame transforms
% This sets Orientations, nOrientations, nSites and AverageOverChi
p_crystalorientations;


% Options parsing and setting.
%---------------------------------------------------------------------

% documented fields
DefaultOptions.Transitions = [];
DefaultOptions.Threshold = 1e-4;
DefaultOptions.Hybrid = 0;
DefaultOptions.HybridNuclei = [];

% undocumented fields
DefaultOptions.nTRKnots = 3;
DefaultOptions.FuzzLevel = 1e-7;
DefaultOptions.MaxKnots = 2000;
DefaultOptions.RediagLimit = 0.95;

DefaultOptions.Intensity = 1;

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

StrainsPresent = any([System.HStrain(:); System.DStrain(:); System.gStrain(:); System.AStrain(:)]);
ComputeStrains = StrainsPresent && (nargout>2);
ComputeIntensities = ((nargout>1) & Opt.Intensity);


% Preparing kernel and perturbing system Hamiltonians.
%-----------------------------------------------------------------------
logmsg(1,'- Preparations');

% :KLUDGE: Add some fuzz to the hyperfine couplings to avoid degeneracies
% if several (equivalent) nuclei are specified.
if (System.nNuclei>1)
  System.A = System.A.*(1 + Opt.FuzzLevel*rand(size(System.A)));
end

CoreSys = System;

% Perturbational treatment of SHF nuclei
if (CoreSys.nNuclei>=1) && Opt.Hybrid
  
  Nucs = nucstringparse(CoreSys.Nucs);
  
  if any(Opt.HybridNuclei>CoreSys.nNuclei)
    error('Opt.HybridNuclei is incorrect!');
  end
  perturbNuclei = ones(1,CoreSys.nNuclei);
  perturbNuclei(Opt.HybridNuclei) = 0;
  
  idx = find(perturbNuclei);
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
        A = System.A(iNuc,idxE);
        R = eye(3);
        if isfield(System,'AFrame')
          R = erot(System.AFrame(iNuc,idxE)).'; % A frame -> molecular frame
        end
        A = R*diag(A)*R.';
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
        Q = [0 0 0];
        R = eye(3);
        if isfield(System,'Q'), Q = System.Q(iNuc,:); end
        if isfield(System,'QFrame')
          R = erot(System.QFrame(iNuc,:)).'; % Q frame -> molecular frame
        end
        Q = R*diag(Q)*R.';
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

if (nPerturbNuclei>0)
  logmsg(1,'  core system with %d spins and %d states',numel(spinvec(CoreSys)),nCore);
  logmsg(1,'  first-order perturbational nuclei with %d states',nSHFNucStates');
else
  if (CoreSys.nNuclei>0)
    logmsg(1,'  full treatment of all nuclei');
  end
end


%=======================================================================
% Construction of transition list
%=======================================================================

UserTransitions = ~isempty(Opt.Transitions);
if (UserTransitions)
  if ischar(Opt.Transitions)
    if strcmp(Opt.Transitions,'all');
      nSStates = prod(2*CoreSys.S+1);
      logmsg(1,'  using all %d transitions',nSStates*(nSStates-1)/2);
      [u,v] = find(triu(ones(nSStates),1));
      Transitions = sortrows([u,v]);
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
  [v,u] = find(tril(ones(nCore),-1));
  Transitions = [u v];
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
logmsg(1,'- Broadenings',nTransitions);
if (ComputeStrains)
  logmsg(1,'  using strains',nTransitions);
  
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
      if any(CoreSys.DFrame(iEl,:))
        R = erot(CoreSys.DFrame(iEl,:)).'; % D frame -> molecular frame
      end
      R = R.'; % molecular frame -> D frame
      
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
      % Compute Hamiltonian derivatives, pre-multiply with strain FWHMs.
      dHdD{iEl} = DeltaF * F_;
      dHdE{iEl} = DeltaG * G_;
      
    end
    clear Sx_ Sy_ Sz_ SxD2_ SyD2_ SzD2_ F_ G_;
  end
  
  % g-A strain
  %-------------------------------------------------
  % g strain tensor is taken to be along the g tensor itself.
  UsegStrain = any(CoreSys.gStrain(:));
  if UsegStrain
    if isfield(CoreSys,'gFrame')
      R = erot(CoreSys.gFrame(1,:)).'; % g frame -> molecular frame
    else
      R = eye(3);
    end
    gStrainMatrix = diag(CoreSys.gStrain./CoreSys.g(1,:));
    gStrainMatrix = R*gStrainMatrix*R.';
  end
  
  UseAStrain = (CoreSys.nNuclei>0) && any(CoreSys.AStrain);
  if UseAStrain
    if isfield(CoreSys,'AFrame')
      R = erot(CoreSys.AFrame(1,:)).'; % A frame -> molecular frame
    else
      R = eye(3);
    end
    
    Ix_ = R(1,1)*sop(CoreSys,2,1)+R(2,1)*sop(CoreSys,2,2)+R(3,1)*sop(CoreSys,2,3);
    Iy_ = R(1,2)*sop(CoreSys,2,1)+R(2,2)*sop(CoreSys,2,2)+R(3,2)*sop(CoreSys,2,3);
    Iz_ = R(1,3)*sop(CoreSys,2,1)+R(2,3)*sop(CoreSys,2,2)+R(3,3)*sop(CoreSys,2,3);
    
    Sx_ = R(1,1)*sop(CoreSys,1,1)+R(1,2)*sop(CoreSys,1,2)+R(1,3)*sop(CoreSys,1,3);
    Sy_ = R(2,1)*sop(CoreSys,1,1)+R(2,2)*sop(CoreSys,1,2)+R(2,3)*sop(CoreSys,1,3);
    Sz_ = R(3,1)*sop(CoreSys,1,1)+R(3,2)*sop(CoreSys,1,2)+R(3,3)*sop(CoreSys,1,3);
    
    dHdAx = CoreSys.AStrain(1)*Ix_*Sx_;
    dHdAy = CoreSys.AStrain(2)*Iy_*Sy_;
    dHdAz = CoreSys.AStrain(3)*Iz_*Sz_;
    
    clear Ix_ Iy_ Iz_ Sx_ Sy_ Sz_
  end
  if any(CoreSys.HStrain), logmsg(2,'  ## using H strain'); end
  if UsegStrain, logmsg(2, ' ## using g strain'); end
  if UseAStrain, logmsg(2,'  ## using A strain'); end
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

startTime = cputime;
logstr = '';

% Calculate transition rates over all orientations (fixed field!).
Pdat = zeros(nTransitions,nOrientations);

if (ComputeIntensities)
  Idat = zeros(nTransitions,nOrientations);
else
  Idat = [];
end

if (ComputeStrains)
  Wdat = zeros(nTransitions,nOrientations);
else
  Wdat = [];
end

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
      if EasySpinLogLevel==1, fprintf(backspace); end
      fprintf(logstr);
    else
      if (nOrientations>1)
        logstr = sprintf('  1/%d orientations, remaining time unknown\n',nOrientations);
        fprintf(logstr);
      end
    end
  end
  
  % Set up Hamiltonians for 3 lab principal directions
  %-----------------------------------------------------
  [xLab,yLab,zLab] = erot(Orientations(iOri,:));
  
  % z laboratoy axis: external static field
  kGzL = zLab(1)*kGxM + zLab(2)*kGyM + zLab(3)*kGzM;
  % x laboratory axis: B1 excitation field
  kGxL = xLab(1)*kGxM + xLab(2)*kGyM + xLab(3)*kGzM;
  % y laboratory vector: needed for integration over all B1 field orientations.
  kGyL = yLab(1)*kGxM + yLab(2)*kGyM + yLab(3)*kGzM;
  
  [Vs,E] = gethamdata(Exp.Field, kF, kGzL);
  Pdat(:,iOri) = E(v) - E(u);
  
  % Calculate intensities if requested
  if (ComputeIntensities)  
        
    for iTrans = nTransitions:-1:1
      Vu = Vs(:,Transitions(iTrans,1)); % lower state (u)
      Vv = Vs(:,Transitions(iTrans,2)); % upper state (Ev>Eu)
      
      % Compute quantum-mechanical transition rate
      mu = [Vv'*kGxL*Vu; Vv'*kGyL*Vu; Vv'*kGzL*Vu];
      if (linearpolarizedMode)
        if (AverageOverChi)
          mu0 = abs(nB0.'*mu);
          TransitionRates(iTrans) = xi1^2*mu0^2 + (1-xi1^2)*(norm(mu)^2-mu0^2)/2;
        else
          TransitionRates(iTrans) = abs(nB1.'*mu)^2;
        end
      elseif (unpolarizedMode)
        if (AverageOverChi)
          mu0 = abs(nB0.'*mu);
          TransitionRates(iTrans) = ((1+xik^2)*norm(mu)^2+(1-3*xik^2)*mu0^2)/4;
        else
          muk = abs(nk.'*mu);
          TransitionRates(iTrans) = (norm(mu)^2-muk^2)/2;
        end
      elseif (circpolarizedMode)
        if (AverageOverChi)
          mu0 = abs(nB0.'*mu);
          imumu = 1i*cross(mu,conj(mu));
          TransitionRates(iTrans) = ((1+xik^2)*norm(mu)^2+(1-3*xik^2)*mu0^2+circpolarizedMode*2*nB0.'*imumu)/4;
        else
          muk = abs(nk.'*mu);
          imumu = 1i*cross(mu,conj(mu));
          TransitionRates(iTrans) = (norm(mu)^2-muk^2+circpolarizedMode*nk.'*imumu)/2;
        end
      end
      
    end
    
    if any(TransitionRates<0)
      logmsg(-inf,'*********** Negative transition rate encountered in resfreqs!! Please report! **********');
    end
    TransitionRates = abs(TransitionRates);
    
    
    % Compute polarizations if temperature is given.
    if (ComputeBoltzmannPopulations)
      
      Populations = ones(nCore,1);
      
      % Pre-factor for thermal equilibrium populations computations.
      BoltzmannPreFactor = -1e6*planck/boltzm/Exp.Temperature; % MHz^-1
      for iState = 1:nCore
        dE = E(iState) - E(1);
        if (dE<1e-10), dE = 0; end % assure we recognize degenerate state even if numerically non-degenerate
        Populations(iState) = exp(BoltzmannPreFactor*dE);
      end
      Populations(isnan(Populations)) = 1;
      Populations = Populations/sum(Populations);
      
      Polarization = Populations(u) - Populations(v);
      if (nPerturbNuclei>0)
        Polarization = Polarization/prod(2*System.I+1);
      end
      Idat(:,iOri) = TransitionRates(:).*Polarization;
      
    else
      % no temperature given
      % same polarization for each electron transition
      %Polarization = Polarization/prod(2*System.S+1); % needed to make consistent with high-temp limit
      Polarization = 1/prod(2*System.I+1);
      Idat(:,iOri) = TransitionRates*Polarization;
    end
  end
  
  % Calculate width if requested.
  %--------------------------------------------------
  if (ComputeStrains)
    for iTrans = 1:nTransitions
      m = @(Op)Vs(:,v(iTrans))'*Op*Vs(:,v(iTrans)) - Vs(:,u(iTrans))'*Op*Vs(:,u(iTrans));
            
      % H strain: Frequency domain residual width tensor.
      %-----------------------------------------------  
      LineWidth2 = (CoreSys.HStrain*zLab.')^2;

      % D strain
      if UseDStrain
        for iEl = 1:CoreSys.nElectrons
          % add D contribution
          LineWidth2 = LineWidth2 + abs(m(dHdD{iEl}))^2;
          % add E contribution
          LineWidth2 = LineWidth2 + abs(m(dHdE{iEl}))^2;
        end
      end
      
      % A strain
      if UseAStrain
        LineWidth2 = LineWidth2 + abs(m(dHdAx))^2;
        LineWidth2 = LineWidth2 + abs(m(dHdAy))^2;
        LineWidth2 = LineWidth2 + abs(m(dHdAz))^2;
      end;
      
      % g strain
      if UsegStrain
        dg2 = (m(kGzL)*zLab*gStrainMatrix*zLab.'*Exp.Field)^2;
        LineWidth2 = LineWidth2 + abs(dg2);
      end
      Wdat(iTrans,iOri) = sqrt(LineWidth2);
    end
  end
  
  % First-order approximation for nuclei
  %-------------------------------------------------------
  if (nPerturbNuclei>0)
    for iTrans = 1 :nTransitions
      U = Vs(:,u(iTrans));
      V = Vs(:,v(iTrans));
      % Calculate <S> for both states involved in the transition
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
if ComputeIntensities && ~UserTransitions
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
    logmsg(2,'  ## all resonances above relative intensity threshold %f',PostSelectionThreshold);
  end
else
  logmsg(2,'  ## no intensities computed, no intensity post-selection');
end

if (EasySpinLogLevel>=2),
  partlyNaN = any(isnan(Pdat),2);
  nChopped = sum(partlyNaN);
  if (nChopped>0)
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
  if (ComputeIntensities), Idat(idxRmv,:) = []; end;
  if (ComputeStrains), Wdat(idxRmv,:) = []; end
  if (nPerturbNuclei>0)
    for iiNuc = 1:nPerturbNuclei
      pPdatN{iiNuc}(idxRmv,:,:) = [];
      pIdatN{iiNuc}(idxRmv,:,:) = [];
    end
  end
end
nTransitions = size(Pdat,1);
logmsg(2,'  ## %2d resonances left',nTransitions);


logmsg(1,'  %d significant transitions with resonances in range',nTransitions);

if (nTransitions==0)
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

% Assert positive intensities and widths
if any(Idat(:)<0),
  logmsg(-inf,'*********** Negative intensity encountered in resfields!! Please report! **********');
end
if any(Wdat(:)<0),
  logmsg(-inf,'*********** Negative width encountered in resfields!! Please report! **************');
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
  if(ComputeStrains)
    if numel(Wdat)>0
      Wdat = reshape(permute(repmat(Wdat,[1,1,nSubTransitions]),[1 3 2]),nTotalTrans,nOrientations);
    end
  end
  Transitions = reshape(permute(repmat(Transitions,[1 1 nSubTransitions]),[1 3 2]),nTotalTrans,2);
  clear pPdat pIdat
end


% Resonance data summary
%---------------------------------------------------------------------
logmsg(2,'  ## resonances min %g MHz, max %g MHz',min(Pdat(:)),max(Pdat(:)));
if (ComputeIntensities)
  logmsg(2,'  ## amplitudes min %g, max %g',min(Idat(:)),max(Idat(:)));
end
if (ComputeStrains && numel(Wdat)>0)
  logmsg(2,'  ## widths min %g mT, max %g mT',min(Wdat(:)),max(Wdat(:)));
end

% Reshape arrays in the case of crystals with site splitting
if (nSites>1) && ~isfield(Opt,'peppercall')
  siz = [nTransitions*nSites, numel(Pdat)/nTransitions/nSites];
  Pdat = reshape(Pdat,siz);
  if ~isempty(Idat), Idat = reshape(Idat,siz); end
  if ~isempty(Wdat), Wdat = reshape(Wdat,siz); end
end

% Arrange the output.
Output = {Pdat,Idat,Wdat,Transitions};
varargout = Output(1:max(nargout,1));

return
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


function [V,E,D1,dE] = gethamdata(B,F,G,idx)

% Compute eigenvalues and eigenvectors of Hamiltonian
[V,E] = eig(F+B*G);
E = diag(E).';

% Compute correct eigenvectors for zero-field degeneracies
if (B==0)
  dE = abs(diff(E)).';
  tol = 1e3*eps*max(dE);
  groupIdx = cumsum([1; dE>tol]);
  GG = V'*G*V;
  GG = (GG+GG')/2; % important: symmetrize
  VV_ = [];
  for iGroup = 1:max(groupIdx)
    ix = find(groupIdx==iGroup);
    [v,e] = eig(GG(ix,ix));
    VV_ = blkdiag(VV_,v);
  end
  V = V*VV_;
end

if (nargout>2)
  % Compute first derivative dE/dB
  D1 = real(diag(V'*G*V)).';
end

if (nargout>3)
  % Compute transition frequencies Ev-Eu
  M = E(:);
  M = M(:,ones(1,length(E)));
  dE = M.' - M;
  dE = dE(idx);
end

return
