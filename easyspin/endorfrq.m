% endorfrq  Compute ENDOR frequencies and intensities 
%
%    Pos = endorfrq(Sys,Exp)
%    Pos = endorfrq(Sys,Exp,Opt)
%    [Pos,Int] = endorfrq(...)
%    [Pos,Int,Tra] = endorfrq(...)
%
%    Sys:  spin system structure
%
%    Exp: experimental parameter structure
%      mwFreq              spectrometer frequency [GHz]
%      Field               magnetic field [mT]
%      Temperature         in K (optional, by default off (NaN))
%      ExciteWidth         FWHM of excitation [MHz] (optional, default Inf)
%      Range               frequency range [MHz] (optional, default [])
%      SampleFrame         Nx3 array of Euler angles (in radians) for sample orientations
%      CrystalSymmetry     crystal symmetry (space group etc.)
%      MolFrame            Euler angles (in radians) for molecular frame orientation
%
%    Opt: computational options structure
%      Verbosity          log level (0, 1 or 2)
%      Transitions        [mx2 array] of level pairs
%      Threshold          for transition selection.
%      Nuclei             index for nuclei (1 to #nuclei)
%      Intensity          'on','off'
%      Enhancement        'on','off'
%      Sites              list of crystal sites to include (default []: all)
%
%    See the documentation for a detailed description of the arguments.

function varargout = endorfrq(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

warning(chkmlver);
if nargin<2 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>4, error('Too many output arguments.'); end

if nargin<3, Opt = struct; end
if isempty(Opt), Opt = struct; end

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

% Initialize optional output structure
Info = struct;

% Process spin system.
%---------------------------------------------------------------------
[Sys,err] = validatespinsys(Sys);
error(err);
if Sys.nNuclei==0
  error('The system doesn''t contain any nuclei.');
end

DefaultSys.lw = 0;
DefaultSys.HStrain = [0 0 0];

Sys = adddefaults(Sys,DefaultSys);

% Process parameters
%---------------------------------------------------------------------
DefaultExp.mwFreq = NaN;
DefaultExp.Temperature = NaN;
DefaultExp.ExciteWidth = Inf;
DefaultExp.Field = NaN;

DefaultExp.SampleFrame = [0 0 0];
DefaultExp.CrystalSymmetry = 1;
DefaultExp.MolFrame = [0 0 0];
DefaultExp.SampleRotation = [];

DefaultExp.Range = []; % for compatibility, internal

if ~isfield(Exp,'Field'), error('Exp.Field is missing.'); end

Exp = adddefaults(Exp,DefaultExp);

mwFreq = Exp.mwFreq*1e3; % GHz -> MHz

computeNonEquiPops = isfield(Sys,'initState') && ~isempty(Sys.initState);
if computeNonEquiPops
  computeBoltzmann = false;
else
  if isinf(Exp.Temperature)
    error('If given, Params.Temperature must be a finite value.');
  end
  computeBoltzmann = ~isnan(Exp.Temperature);
end

if isfield(Exp,'ExciteWidth')
  if ~isfield(Exp,'mwFreq'), error('Par.mwFreq is missing.'); end
end

if ~isfield(Opt,'Sites')
  Opt.Sites = [];
end

% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,nSites,AverageOverChi] = p_crystalorientations(Exp,Opt);

%-----------------------------------------------------------------------


% Process Options.
%-----------------------------------------------------------------------
% Obsolete fields
ObsoleteOptions = {'nTransitions'};
for k = 1:numel(ObsoleteOptions)
  if isfield(Opt,ObsoleteOptions{k})
    error('Options.%s is obsolete. Please remove from code!',ObsoleteOptions{k});
  end
end

% Documented fields
DefaultOpt.Transitions = [];
% Opt.Threshold < 1e-3 sometimes results in a forbidden EPR
% transition being included in the Transitions list. But it
% should not harm, since its frequency is much higher and can
% easily be removed.
DefaultOpt.Threshold = 1e-3;
DefaultOpt.Nuclei = [];
DefaultOpt.Intensity = 'on';
% Opt.Enhancement is off by default, since in many ENDOR experiments
% RF transmitter characteristics cancel the effect.
DefaultOpt.Enhancement = 'off';

% Undocumented fields
DefaultOpt.OriWeights = [];
DefaultOpt.OriThreshold = 1e-4;
DefaultOpt.TPSGridSize = 4;      % grid size for transition pre-selection
DefaultOpt.TPSGridSymm = 'D2h';  % grid symmetry for transition pre-selection

%DefaultOpt.PreSelection = 1;
DefaultOpt.SelectionThreshold = 0.01;

Opt = adddefaults(Opt,DefaultOpt);

% List of nuclei whose Zeeman interaction should be included
% for the detection operator in the transition selection.
if isempty(Opt.Nuclei), Opt.Nuclei = 1:Sys.nNuclei; end
if any(Opt.Nuclei<1) || any(Opt.Nuclei>Sys.nNuclei)
  error('Opt.Nuclei is out of range.');
end

IntensitySwitch = parseoption(Opt,'Intensity',{'off','on'}) - 1;

EnhancementSwitch = parseoption(Opt,'Enhancement',{'off','on'}) - 1;

UseOriWeights = ~isempty(Opt.OriWeights) & (Opt.OriThreshold>0);
if UseOriWeights
  if numel(Opt.OriWeights)~=nOrientations
    error('Opt.OriWeights has wrong number of elements (%d).\n It must have one entry for each of the %d orientations!',...
    numel(Opt.OriWeights),nOrientations);
  end
end

% Preparing Hamiltonian representations etc.
%-----------------------------------------------------------------------
% The first and only compilation of the full Hamiltonian.
[H0,muxM,muyM,muzM] = ham(Sys);


% For polarized systems, pre-compute ZF eigenstates.
if computeNonEquiPops

  initState = Sys.initState{1};
  initStateBasis = Sys.initState{2};

  % Check and adapt input dimensions
  nElectronStates = prod(2*Sys.S+1);
  nNucStates = prod(Sys.I*2+1);
  nStates = nElectronStates*nNucStates;

  [sz1,sz2] = size(initState);
  if sz1==sz2
    % Density matrix
    if nElectronStates~=sz1 && nStates~=sz1
      error('The density matrix in Sys.initState must have dimensions of nxn with n = %d or %d.',nElectronStates,nStates)
    end
    if numel(initState)==nElectronStates^2 && numel(initState)~=nStates^2
      initState = kron(initState,eye(nStates/nElectronStates));
    end
    initState = initState/trace(initState);
  else
    % Vector of populations
    if numel(initState)~=nElectronStates && numel(initState)~=nStates
      error('The population vector in Sys.initState must have %d or %d elements.',nElectronStates,nStates);
    end
    initState = initState(:);
    if numel(initState)==nElectronStates && numel(initState)~=nStates
      initState = kron(initState,ones(nStates/nElectronStates,1));
    end
    initState = initState/sum(initState);
    % Convert population vector to density matrix for populations provided in eigenbasis
    if strcmp(initStateBasis,'eigen')
      initState = diag(initState);
    end
  end

  if strcmp(initStateBasis,'zerofield')
    % Pre-compute zero-field energies and eigenstates
    [ZFStates,ZFEnergies] = eig(H0);
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
  end

else
%   ZFEnergies = sort(real(eig(F)));
end

% Intensities are computed if option settings are positive and
% if the user has provided an output variable.
ComputeIntensities = (nargout>1) & (IntensitySwitch);

% Specifies whether to compute eigenvectors. The variable is
% redundant at the moment, since only intensity computations
% need eigenvectors.
ComputeVectors = ComputeIntensities;

OrientationSelection = isfinite(Exp.ExciteWidth) & (Exp.ExciteWidth>0);

if OrientationSelection
  logmsg(2,'  including excitation width %g MHz',Exp.ExciteWidth);
end

nStates = Sys.nStates; % state space dimension

%=======================================================================
%                         Transition pre-selection
%=======================================================================
% The following part composes a list of transitions, ie level pairs 
% for which peak data are computed. The list either is chosen from
% a user specified list in Opt.Transitions, or is composed by an
% automatic procedure which selects the most intense transitions,
% their number being determined by a threshold for the relative transitions
% rate Opt.Threshold. The relative transition rate for the most intense
% transition is equal to 1.

if isempty(Opt.Transitions)
  logmsg(1,'  automatic transition selection');
  logmsg(2,'    (threshold %g, grid size %d, grid symmetry %s)',Opt.Threshold,Opt.TPSGridSize,Opt.TPSGridSymm);
  
  % Set a coarse grid, independent of the effective symmetry of
  % the Hamiltonian.
  grid = sphgrid(Opt.TPSGridSymm,Opt.TPSGridSize);
  phi = grid.phi;
  theta = grid.theta;
  TRWeights = grid.weights;
  
  % pre-compute trigonometric functions
  st = sin(theta);
  sp = sin(phi);
  ct = cos(theta);
  cp = cos(phi);
  
  % Compute selection detection operators (NMR transitions only!)
  [smuxM,smuyM,smuzM] = ham_nz(Sys,Opt.Nuclei);
  
  % preallocate the transition rate matrix
  TransitionRates = zeros(nStates);
  maxE = zeros(nStates);
  minE = Inf(nStates);
  
  % calculate transition rates over all orientations
  for iOri = 1:length(theta)
    
    % solve eigenproblem
    [Vs,E] = eig(H0 - Exp.Field*...
      (st(iOri)*(cp(iOri)*muxM + sp(iOri)*muyM) + ct(iOri)*muzM));
    E = diag(E);
    
    % sum up transition rates
    smupM = cp(iOri)*smuxM+ sp(iOri)*smuyM;
    TransitionRates = TransitionRates + TRWeights(iOri) * ...
      (abs(Vs'*( ct(iOri)*smupM - st(iOri)*smuzM)*Vs).^2 + ...
       abs(Vs'*(-sp(iOri)*smuxM + cp(iOri)*smuyM)*Vs).^2)/2;
       
    % Determine minima and maxima of transition frequencies
    EE = E(:,ones(1,nStates));
    EE = abs(EE - EE.');
    idx = EE > maxE; maxE(idx) = EE(idx);
    idx = EE < minE; minE(idx) = EE(idx);
  end
  
  % Remove unused matrices. With high nStates, they take a lot of space.
  clear Vs E EE smupM smuxM smuyM smuzM st ct sp cp TRWeights idx;
  
  % Remove transitions completely out of range.
  if ~isempty(Exp.Range) && ~any(isnan(Exp.Range))
    TransitionRates((maxE<Exp.Range(1))|(minE>Exp.Range(2))) = 0;
  end
  clear maxE minE;
  
  % Remove lower triangular part and electronic transitions
  nSStates = prod(2*Sys.S+1);
  idxElectronicTransitions = ~kron(eye(nSStates),ones(nStates/nSStates));
  idxLowerTriangle = logical(tril(ones(nStates)));
  keepidx = ~(idxLowerTriangle | idxElectronicTransitions);
  TransitionRates(~keepidx) = [];
  clear idxLowerTriangle idxElectronicTransitions;
  
  [u,v] = find(keepidx); % Compute level pair indices.
  Transitions = [u,v];
  clear keepidx u v;

  % Use threshold for number determination.
  nTransitions = sum(TransitionRates>Opt.Threshold*max(TransitionRates));
  
  % Sort TransitionRates in descending order!
  [~,idx] = sort(-TransitionRates);
  % Select most intense transitions.
  Transitions = Transitions(idx(1:nTransitions),:);
  clear idx unused TransitionRates;
  
else % User-specified transitions present.
  logmsg(1,'  using %d user-specified transitions',size(Opt.Transitions,1));
  % Guarantee that lower index comes first (gives later u < v).
  Transitions = sort(Opt.Transitions,2);
end

u = Transitions(:,1);
v = Transitions(:,2);
nTransitions = length(u);
TRidx = u + (v-1)*nStates; % Indices into UPPER triangle.

% Diagnostic display.
logmsg(1,'  %d transitions selected',nTransitions);
if EasySpinLogLevel>=3, disp(Transitions); end

% Issue a warning if the resulting transition list is empty.
if isempty(Transitions)
  warning('No transitions found. Decrease Opt.Threshold.');
end



%=======================================================================
%                      PEAK DATA GENERATION
%=======================================================================

% Pre-allocations and initializations
%-----------------------------------------------------------------------
% Display information on what is going to be computed.
msg = '  computing peak positions';
if (ComputeIntensities), msg = [msg ', intensities']; end
logmsg(1,msg);

% Preallocations.
Pdat = NaN(nTransitions,nOrientations);
Idat = [];
if ComputeIntensities, Idat = zeros(nTransitions,nOrientations); end

% Other preparations.
if ComputeIntensities
  % Set detection operators for intensity computations.
  [DxM,DyM,DzM] = ham_nz(Sys,Opt.Nuclei);
  if EnhancementSwitch
    % Zeeman interaction including electronic Zeeman interaction,
    % includes implicitely hyperfine enhancement
    [DxMe,DyMe,DzMe] = ham_ez(Sys);
    DxM = DxM + DxMe;
    DyM = DyM + DyMe;
    DzM = DzM + DzMe;
  end
  
  if OrientationSelection
    % Electron Zeeman interaction operators for EPR transition
    % rate computation
    [ExM,EyM,EzM] = ham_ez(Sys,1);
    % pre-square line width
    lwExcite2 = Exp.ExciteWidth^2;
  end
  
  % Prefactor used for computing the Boltzmann population distribution.
  if computeBoltzmann
    BoltzmannPreFactor = -1e6*planck/boltzm/Exp.Temperature;
  end
end

% Initialize parameters for orientation selectivity determination.
% Selectivity = (maxEPRfreq-minEPRfreq)/minExWidth
if OrientationSelection
  maxEPRfreq = -Inf;
  minEPRfreq = Inf;
  minExWidth = Inf;
end

% Keep only orientations with weights above weight threshold.
if UseOriWeights
  logmsg(1,'  user-supplied orientation pre-selection: skipping %d of %d orientations',...
    sum(Opt.OriWeights<Opt.OriThreshold),nOrientations);
end

% Loop over all (remaining) orientations.
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
      fprintf([backspace logstr]);
    else
      if nOrientations>1
        logstr = sprintf('  1/%d orientations, remaining time unknown\n',nOrientations);
        fprintf(logstr);
      end
    end
  end
  
  if UseOriWeights
    if Opt.OriWeights(iOri)<Opt.OriThreshold; continue; end
  end
  
  % Lab frame axes in molecular frame representation.
  [xLab,yLab,zLab] = erot(Orientations(iOri,:),'rows');
  
  
  % Compute eigenstate energies and vectors.
  H = H0 - Exp.Field*(zLab(1)*muxM + zLab(2)*muyM + zLab(3)*muzM);
  if ComputeVectors
    [Vs,E0] = eig(H);
    E0 = sort(diag(E0));
  else
    E0 = sort(eig(H));
  end
  
  % Position data
  Pdat(:,iOri) = E0(v) - E0(u);
  
  if ComputeIntensities
      
    % mw: xLab, RF: yLab
    DyL = yLab(1)*DxM + yLab(2)*DyM + yLab(3)*DzM;
    if (AverageOverChi)
      % average over lab xy plane
      DxL = xLab(1)*DxM + xLab(2)*DyM + xLab(3)*DzM;
      EndorIntensity = (abs(Vs'*DxL*Vs).^2 + abs(Vs'*DyL*Vs).^2)/2;
      if OrientationSelection
        ExL = xLab(1)*ExM + xLab(2)*EyM + xLab(3)*EzM;
        EyL = yLab(1)*ExM + yLab(2)*EyM + yLab(3)*EzM;
        EPRIntensity = (2*pi)*(abs(Vs'*ExL*Vs).^2+abs(Vs'*EyL*Vs).^2)/2;
      end
    else
      EndorIntensity = abs(Vs'*DyL*Vs).^2;
      if OrientationSelection
        ExL = xLab(1)*ExM + xLab(2)*EyM + xLab(3)*EzM;
        EPRIntensity = abs(Vs'*ExL*Vs).^2;
      end
    end
    
    % Pulse ENDOR: include rf flip angle
    if ~isfield(Exp,'rf'), Exp.rf = 0; end
    pulseENDOR = Exp.rf>0;
    if pulseENDOR
      EndorIntensity = (1 - cos(sqrt(EndorIntensity)*Exp.rf))/2;
    end
    
    % Compute excitation factor
    if OrientationSelection
      lw2 = sum((Sys.HStrain(:).*zLab).^2) + lwExcite2;
      
      % weighting with EPR line width and excitation band width
      dE = E0(:,ones(1,nStates));
      dE = dE.' - dE; % positive in upper triangle
      
      % Compute polarization if temperature or zero-field populations are given.
      Populations = [];
      if computeBoltzmann
        Populations = exp(BoltzmannPreFactor*(E0-E0(1)));
        %Polarization = (Populations(u) - Populations(v))/sum(Populations);
      elseif computeNonEquiPops
        %Populations = (abs(ZFStates'*Vs).^2).'*ZFPopulations; % lower level
        %Polarization = PopulationU - PopulationV;
        Populations = zeros(size(dE));
        switch initStateBasis
          case 'eigen'
            for iState = 1:nStates
              Populations(iState) = initState(iState,iState);
            end
          otherwise
            for iState = 1:nStates
              Populations(iState) = Vs(:,iState)'*initState*Vs(:,iState);
            end
        end
      end
      if isempty(Populations)
        Polarizations = ones(size(dE));
      else
        P = Populations(:,ones(1,nStates));
        Polarizations = abs(P - P.');
      end
      
      % weight EPR transition probabilities with excitation shape
      EPRIntensity = sqrt(2/pi)*...
        Polarizations.*EPRIntensity.*exp(-2/lw2*(mwFreq-abs(dE)).^2);

      % sum up tp for all transitions leading to a specific level,
      % exclude the level itself
      SumTP = (sum(EPRIntensity) - diag(EPRIntensity).').';
      EPRIntensity = SumTP(u) + SumTP(v) - 2*EPRIntensity(TRidx);
      
      % select only needed intensities
      ExcitationFactor = EPRIntensity;
      
      % orientation selectivity
      dE(dE<0.3*mwFreq) = [];
      minFreq = min(dE);
      maxFreq = max(dE);
      if minFreq<minEPRfreq, minEPRfreq = minFreq; end
      if maxFreq>maxEPRfreq, maxEPRfreq = maxFreq; end
      if lw2<minExWidth^2, minExWidth = sqrt(lw2); end

    else
      ExcitationFactor = 1;
    end % if OrientationSelection else
        
    % Compute polarization if temperature or zero-field populations are given.
    if computeBoltzmann
      Populations = exp(BoltzmannPreFactor*(E0-E0(1)));
      NuclearPolarization = (Populations(u) - Populations(v))/sum(Populations);
    else
      NuclearPolarization = 1;
    end

    % Total intensity
    Idat(:,iOri) = NuclearPolarization .* ExcitationFactor .* EndorIntensity(TRidx);  %#ok
    
  end % if ComputeIntensity
  
end % for all orientations
%=======================================================================

% Compute selectivity.
if (OrientationSelection)
  Selectivity = (maxEPRfreq-minEPRfreq)/minExWidth;
  logmsg(2,'  limited excitation width, selectivity %g',Selectivity);
else
  Selectivity = 0;
  logmsg(2,'  infinite excitation width');
end

% Information structure
Info.Selectivity = Selectivity;

% Reshape arrays in the case of crystals with site splitting
d = dbstack;
saltCall = numel(d)>2 && strcmp(d(2).name,'salt');
if nSites>1 && ~saltCall
  siz = [nTransitions*nSites, numel(Pdat)/nTransitions/nSites];
  Pdat = reshape(Pdat,siz);
  if ~isempty(Idat), Idat = reshape(Idat,siz); end
end

% Arrange output.
Output = {Pdat,Idat,Transitions,Info};
varargout = Output(1:nargout);
if nargout==0, varargout = {Pdat};end
end
