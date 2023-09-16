% saffron    Simulate pulse EPR signals and spectra
%
%     S = saffron(Sys,Exp,Opt)
%     [x,S] = saffron(Sys,Exp,Opt)
%     [x,S,info] = saffron(Sys,Exp,Opt)
%
%  Inputs:
%     Sys   ... spin system with electron spin and ESEEM nuclei
%     Exp   ... experimental parameters (time unit µs)
%     Opt   ... simulation options
%
%  Outputs:
%     x       ... axis/axes for S, contains all indirect dimensions and for
%                 transient detection transient time axis of transient, for
%                 ENDOR the frequency axis
%     S       ... simulated signal (ESEEM) or spectrum (ENDOR)
%     info    ... structure with FFT of ESEEM signal etc.

function varargout = saffron(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version
error(chkmlver);

% Get time for performance report at the end
startTime = datetime;

% Input argument scanning, get display level and prompt
%===============================================================================

% Guard against wrong number of input or output arguments
if nargin<1, error('Please supply a spin system as first input argument.'); end
if nargin<2, error('Please supply experimental parameters as second input argument.'); end
if nargin>3, error('Too many input arguments, the maximum is three.'); end
if nargout>3, error('Too many output arguments.'); end

% Initialize options structure if not given
if nargin<3, Opt = struct; end
if isempty(Opt), Opt = struct; end

if ~isstruct(Sys) && ~iscell(Sys)
  error('The first input (Sys) must be a structure or a cell array of structures.');
end
if ~isstruct(Exp)
  error('The second input argument (Exp) must be a structure.');
end
if ~isstruct(Opt)
  error('The third input argument (Opt) must be a structure.');
end

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

% Process Opt.separate
if ~isfield(Opt,'separate'), Opt.separate = ''; end
[separateOutput,err] = parseoption(Opt,'separate',{'','components','transitions','orientations','sites'});
error(err);
separateComponentSpectra = separateOutput==2;
separateTransitionSpectra = separateOutput==3;
separateOrientationSpectra = separateOutput==4;
separateSiteSpectra = separateOutput==5;
if separateTransitionSpectra || separateSiteSpectra || separateOrientationSpectra
  separateComponentSpectra = true;
end
if separateTransitionSpectra
  error('saffron does not support Opt.separate=''transitions''.')
end
if separateSiteSpectra
  error('saffron does not support Opt.separate=''sites''.')
end
if separateOrientationSpectra
  error('saffron does not support Opt.separate=''orientations''.')
end

%===============================================================================
% Loop over species and isotopologues
%===============================================================================
isENDOR = ischar(Exp.Sequence) && Exp.Sequence=="MimsENDOR";

if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end

singleIsotopologue = isfield(Sys,'singleiso') && Sys.singleiso;
if ~singleIsotopologue

  thirdOutput = nargout>=3 | nargout==0;
  autoRange = false;
  [x,data,info] = compisoloop(@saffron,Sys,Exp,Opt,autoRange,thirdOutput,separateComponentSpectra);

  switch nargout
    case 0
      if isfield(info,'td')
        saffron_plot(x,info,isENDOR,Sys,Exp,Opt);
      else
        Opt.SinglePointDetection = info.SinglePointDetection;
        s_plotting(x,data,Exp,Opt);
      end
    case 1
      varargout = {data};
    case 2
      varargout = {x,data};
    case 3
      varargout = {x,data,info};
  end
  return
end


%===============================================================================
% Single-isotopologue simulation
%===============================================================================

logmsg(1,['=begin=saffron====' char(datetime) '=================']);
logmsg(2,'  log level %d',EasySpinLogLevel);
logmsg(1,'-general-----------------------------------------------');


%===============================================================================
% Spin system structure
%===============================================================================
[Sys,err] = validatespinsys(Sys);
error(err);
if Sys.MO_present
  error('saffron does not support Sys.Ham* parameters.');
end
if any(Sys.L(:))
  error('saffron does not support Sys.L.');
end

% Error on spidyan-specific fields
spidyanSysFields = {'initState','eqState'};
for f = 1:numel(spidyanSysFields)
  field_ = spidyanSysFields{f};
  if isfield(Sys,field_) && ~isempty(Sys.(field_))
    error('Sys.%s is specific to spidyan and cannot be used with saffron. Please remove.',...
      field_);
  end
end


%===============================================================================
% Experiment structure
%===============================================================================

DefaultExp.Temperature = [];
DefaultExp.Ordering = [];
DefaultExp.SampleFrame = [];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];
DefaultExp.SampleRotation = [];
Exp = adddefaults(Exp,DefaultExp);

% Check for obsolete fields and unsupported fields
if isfield(Exp,'Orientations')
  error('Exp.Orientations is no longer supported, use Exp.SampleFrame instead.');
end
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame instead.');
end
if isfield(Exp,'lightBeam') && ~isempty(Exp.lightBeam)
  error('Photoselection (via Exp.lightBeam) is not supported.')
end

% Field
if ~isfield(Exp,'Field')
  error('Exp.Field is missing. Give a magnetic field in mT.');
end

% Exp.Sequence is a required field
if ~isfield(Exp,'Sequence')
  error('Exp.Sequence is missing. Provide a pre-defined or user-defined pulse sequence.');
end

% Temperature
if ~isempty(Exp.Temperature)
  error('Exp.Temperature is not supported for pulse EPR simulations.');
end

% Powder vs. crystal simulation
PowderSimulation = isempty(Exp.MolFrame) && isempty(Exp.CrystalSymmetry);
Exp.PowderSimulation = PowderSimulation;

% Partial ordering
if ~isempty(Exp.Ordering)
  if ~PowderSimulation
    error('Partial ordering (Exp.Ordering) can only be used in a powder simulation.');
  else
    error('Partial ordering (Exp.Ordering) is not implemented in saffron.');
  end
end


if isfield(Exp,'nPoints')
  nDataPoints = prod(Exp.nPoints);
else
  nDataPoints = 1;
end

%===============================================================================
% Options structure
%===============================================================================
if ~isfield(Opt,'GridSymmetry'), Opt.GridSymmetry = []; end
if ~isfield(Opt,'GridFrame'), Opt.GridFrame = []; end
if ~isfield(Opt,'Transitions'), Opt.Transitions = []; end
if ~isfield(Opt,'Sites'), Opt.Sites = []; end

% Obsolete options
if isfield(Opt,'nKnots')
  error('Options.nKnots is obsolete. Use Options.GridSize instead, e.g. Options.GridSize = 91.');
end

% Unsupported/forbidden options
spidyanOptFields = {'StateTrajectories','ExcOperator'};
for index = 1 : numel(spidyanOptFields)
  field_ = spidyanOptFields{index};
  if isfield(Opt,field_)
    Opt = rmfield(Opt,field_);
    warning(['Opt.' field_ ' is specific to spidyan and can not be used with saffron. It will be ignored.'])
  end
end

% Relaxation
if isfield(Opt,'Relaxation') && length(Opt.Relaxation) > 1
  error('Opt.Relaxation has to be true or false.')
elseif ~isfield(Opt,'Relaxation') && (isfield(Sys,'T1') || isfield(Sys,'T2'))
  Opt.Relaxation = true;
elseif ~isfield(Opt,'Relaxation')
  Opt.Relaxation = false;
end

% Grid size and interpolation
if ~isfield(Opt,'GridSize'), Opt.GridSize = 30+1; end
if numel(Opt.GridSize)>1
  error('Only one number allowed in Opt.GridSize. saffron does not support interpolation.');
end
if Opt.GridSize<7
  error('Opt.GridSize must be at least 7. You gave %d.',Opt.GridSize);
end

% Determine simulation mode
if ~isfield(Opt,'SimulationMode'), Opt.SimulationMode = 'fast'; end

% Parse pulse sequence
[Events,Vary,Opt,oldSeq] = s_sequencer(Exp,Opt);
if ~isempty(oldSeq)
  Exp.Flip = oldSeq.Flip;
  Exp.Phase = oldSeq.Phase;
  Exp.tp = oldSeq.tp;
  Exp.t = oldSeq.t;
  Exp.Inc = oldSeq.Inc;
  Exp.dt = oldSeq.dt;
end

[~,err] = parseoption(Opt,'SimulationMode',{'fast','thyme'});
error(err);
fastSimulationMode = Opt.SimulationMode=="fast";

if ~fastSimulationMode  && ~nargout && isfield(Exp,'nPoints')
  if Opt.SinglePointDetection  && length(Exp.nPoints) > 2
    error('Not enough output arguments for more than two indirect dimensions with single point detections. The data can not be generally plotted.')
  elseif ~Opt.SinglePointDetection  && length(Exp.nPoints) > 1
    error('Not enough output arguments for more than one indirect dimensions with transient detections. The data can not be generally plotted.')
  end
end

if fastSimulationMode
  if Sys.nElectrons>1
    error('Fast mode does not work for more than one electron spin.');
  end
  if isfield(Sys,'n') && any(Sys.n~=1)
    error('Fast mode does not work for equivalent nuclei.');
  end
  if isfield(Sys,'nn') && any(Sys.nn(:)~=0)
    error('Fast mode does not work in the presence of nuclear-nuclear couplings.');
  end
end


% Set up orientation loop
Exp.R_sample = p_samplerotmatrix(Exp.SampleRotation);
if ~isfield(Exp,'OriWeights')
  [Exp,Opt] = p_symandgrid(Sys,Exp,Opt);
end
% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,nSites] = p_crystalorientations(Exp,Opt);
if numel(Exp.OriWeights)~=nOrientations
  Exp.OriWeights = repmat(Exp.OriWeights,1,nSites);
end



if fastSimulationMode

  %===================================================================
  % Experiment structure
  %===================================================================

  spidyanExpFields = {'DetSequence','DetOperator'};
  if any(isfield(Exp,spidyanExpFields))
    error('Exp.DetSequence and Exp.DetOperator are specific to spidyan and can not be used with saffron. Please remove them.');
  end

  % T1, T2
  if isfield(Sys,'T1T2')
    error('T1 and T2 should be provided in Sys.T1 and Sys.T2 and not in Sys.T1T2.');
  end
  Sys.T1T2 = [0 0];
  if isfield(Sys,'T1')
    if any(size(Sys.T1) > 1)
      error('The fast algorithm does not support transition selective relaxation times, please provide a single value for Sys.T1.')
    else
      Sys.T1T2(1) = Sys.T1;
    end
  end
  if isfield(Sys,'T2')
    if any(size(Sys.T2) > 1)
      error('The fast algorithm does not support transition selective relaxation times, please provide a single value for Sys.T2.')
    else
      Sys.T1T2(2) = Sys.T2;
    end
  end
  Sys.T1T2(Sys.T1T2==0) = inf;

  if any(Sys.T1T2<=0) || any(~isreal(Sys.T1T2))
    error('T1 and T2 in Sys.T1 and Sys.T2 must be positive, in microseconds.');
  end

  % Pulse sequence
  predefinedExperiment = isfield(Exp,'Sequence') && ~isempty(Exp.Sequence) && ischar(Exp.Sequence);
  if predefinedExperiment

    if isfield(Exp,'Filter')
      error('Exp.Filter can only be used with custom sequences.');
    end

    ExperimentNames = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR'};
    Exp.ExperimentID = find(strcmp(Exp.Sequence,ExperimentNames));
    if isempty(Exp.ExperimentID)
      error('Exp.Sequence ''%s'' not recognized.',Exp.Sequence);
    end
    if numel(Exp.ExperimentID)>1, error('Ambiguous sequence name.'); end
    logmsg(1,'Sequence: %s',ExperimentNames{Exp.ExperimentID});

    if isfield(Exp,'tp')
      if any(Exp.tp~=0)
        error('You cannot use predefined sequences (Exp.Sequence) with real pulses (Exp.tp).');
      end
    end

    isENDOR = false;
    tauRequired = true;
    switch Exp.ExperimentID
      case 1  % 2pESEEM
        nIntervals = 2; nDimensions = 1; IncSchemeID = 2;
        nPathways = 1; pathwayprefactor = +1/2;
        tauRequired = false;
      case 2  % 3pESEEM
        nIntervals = 3; nDimensions = 1; IncSchemeID = 1;
        nPathways = 2; pathwayprefactor = +1/8*[1 1];
      case 3  % 4pESEEM
        nIntervals = 4; nDimensions = 1; IncSchemeID = 2;
        nPathways = 2; pathwayprefactor = -1/8*[1 1];
      case 4  % HYSCORE
        nIntervals = 4; nDimensions = 2; IncSchemeID = 11;
        nPathways = 2; pathwayprefactor = -1/8*[1 1];
      case 5  % Mims ENDOR
        isENDOR = true;
        nIntervals = 3; nDimensions = 1; IncSchemeID = 0;
        nPathways = 2; pathwayprefactor = +1/8*[1 1];
    end

    if tauRequired
      if ~isfield(Exp,'tau')
        error('Exp.tau is missing.');
      end
      if numel(Exp.tau)~=1
        error('Exp.tau must contain a single value (in µs).');
      end
      if Exp.tau==0
        error('Exp.tau must be larger than 0.');
      end
    else
      if ~isfield(Exp,'tau')
        Exp.tau = 0;
      end
      if numel(Exp.tau)~=1
        error('Exp.tau must contain a single value (in µs).');
      end
    end

  else

    % User-specified pulse sequence -----------------------------------------
    logmsg(1,'User-specified pulse experiment.');
    Exp.ExperimentID = -1;
    isENDOR = false;

    if any(~isinf(Sys.T1T2))
      error('Relaxation times T1 and T2 for custom sequences not supported.');
    end

    if ~isfield(Exp,'Flip')
      error('Fast sequence engine: the pulse flip angles (Exp.Flip) are missing.');
    end
    nIntervals = numel(Exp.Flip);

    % Incrementation scheme
    if ~isfield(Exp,'Inc')
      error('Fast sequence engine: the incrementation scheme (Exp.Inc) is missing.');
    end
    if numel(Exp.Inc)~=nIntervals
      error('Exp.Inc must contain the same number of elements as Exp.Flip.');
    end
    IncScheme = Exp.Inc(Exp.Inc~=0);
    if isequal(IncScheme,1), IncSchemeID = 1;
    elseif isequal(IncScheme,[1 1]), IncSchemeID = 2;
    elseif isequal(IncScheme,[1 -1]), IncSchemeID = 3;
    elseif isequal(IncScheme,[1 2]), IncSchemeID = 11;
    elseif isequal(IncScheme,[1 2 1]), IncSchemeID = 12;
    elseif isequal(IncScheme,[1 2 2]), IncSchemeID = 13;
    elseif isequal(IncScheme,[1 1 2]), IncSchemeID = 14;
    elseif isequal(IncScheme,[1 2 2 1]), IncSchemeID = 15;
    elseif isequal(IncScheme,[1 2 -2 1]), IncSchemeID = 16;
    elseif isequal(IncScheme,[1 1 2 2]), IncSchemeID = 17;
    else
      error('Unsupported incrementation scheme in Exp.Inc.');
    end
    nDimensions = max(abs(Exp.Inc));

    if ~isfield(Exp,'t')
      Exp.t = zeros(1,nIntervals);
    end
    if numel(Exp.t)~=nIntervals
      error('Exp.t must contain the same number of elements as Exp.Flip');
    end
    if all(Exp.t==0) && (numel(IncScheme)<nIntervals)
      logmsg(0,'Some delays are zero, but are not incremented!');
    end

  end

  % dt
  if ~isENDOR
    if ~isfield(Exp,'dt')
      error('Exp.dt is missing.');
    end
    if numel(Exp.dt)==1
      Exp.dt = Exp.dt*ones(1,nDimensions);
    elseif numel(Exp.dt)~=nDimensions
      error('Exp.dt needs either 1 or %d elements, one per dimension. You gave %d.',nDimensions,numel(Exp.dt));
    end
  end

  % nPoints
  if ~isfield(Exp,'nPoints')
    if isENDOR
      Exp.nPoints = 1001;
    else
      if nDimensions==1
        Exp.nPoints = 512;
      else
        Exp.nPoints = [1 1]*256;
      end
    end
  end
  if numel(Exp.nPoints)==1
    Exp.nPoints = Exp.nPoints*ones(1,nDimensions);
  elseif numel(Exp.nPoints)~=nDimensions
    error('Exp.nPoints needs either 1 or %d elements, one per dimension. You gave %d.',nDimensions,numel(Exp.nPoints));
  end

  % tp
  % Determine whether to use ideal pulse theory
  if ~isfield(Exp,'tp')
    Exp.tp = zeros(1,nIntervals);
  else
    if numel(Exp.tp)~=nIntervals
      error('Exp.tp contains a wrong number of elements.');
    end
  end
  idealPulse = (Exp.tp==0);
  realPulse = ~idealPulse;

  % excitation width
  if any(realPulse)
    if isfield(Exp,'ExciteWidth')
      %error('Cannot Exp.ExciteWidth and real pulses (Exp.tp) at the same time.')
    end
  end


  % User-defined experiment: pathways etc.
  %--------------------------------------------------------------------------

  if ~predefinedExperiment

    % Validate incrementation scheme: assert all delays are nonnegative
    for iInterval = 1:nIntervals
      iDim = abs(Exp.Inc(iInterval));
      if iDim==0  % delay is kept constant
        t_range = Exp.t(iInterval)*[1 1];
      else % delay is incremented/decremented
        t_range = Exp.t(iInterval) + sign(Exp.Inc(iInterval))*(Exp.nPoints(iDim)-1)*Exp.dt(iDim);
      end
      if any(t_range<0)
        error('Negative delay after pulse %d.',iInterval);
      end
    end

    % Determine pathways contributing to the echo
    logmsg(1,'  determining pathways contributing to the echo');
    if isfield(Exp,'Pathways')
      code0('ab+-') = [1 2 3 4];
      pathwayList = code0(Exp.Pathways);
      if size(pathwayList,2)~=nIntervals
        error('Exp.Pathways contains pathways longer or shorter than the pulse sequence.');
      end
    else
      pathwayList = sf_pathways(Exp.t,Exp.Inc);
    end
    nPathways = size(pathwayList,1);
    if nPathways==0
      error('Sorry, no focused echo with this sequence and timings.');
    end

    % Apply coherence filters ---------------------
    if ~isfield(Exp,'Filter')
      Exp.Filter = [];
    end
    if ~isempty(Exp.Filter)
      logmsg(1,'  applying user-supplied coherence filters');
      if ~ischar(Exp.Filter)
        error('Exp.Filter must be a string containing ''0'', ''1'', ''a'', ''b'', ''+'', ''-'' and/or ''.''.');
      end
      if numel(Exp.Filter)~=nIntervals
        error('Exp.Filter must have %d elements instead of the given %d elements.',nIntervals,numel(Exp.Filter));
      end
      for iInt = 1:nIntervals
        switch Exp.Filter(iInt)
          case '0', keep = (pathwayList(:,iInt)==1) | (pathwayList(:,iInt)==2);
          case '1', keep = (pathwayList(:,iInt)==3) | (pathwayList(:,iInt)==4);
          case 'a', keep = (pathwayList(:,iInt)==1);
          case 'b', keep = (pathwayList(:,iInt)==2);
          case '+', keep = (pathwayList(:,iInt)==3);
          case '-', keep = (pathwayList(:,iInt)==4);
          case '.', keep = (pathwayList(:,iInt)~=0);
          case '*', keep = (pathwayList(:,iInt)~=0);
          otherwise
            error('Exp.Filter(%d) = ''%s'' is invalid. Use ''0'', ''1'', ''a'', ''b'', ''+'', ''-'' or ''.''.',iInt,Exp.Filter(iInt));
        end
        pathwayList = pathwayList(keep,:);
      end
      nPathways = size(pathwayList,1);
      if nPathways==0
        error('Exp.Filter is too restrictive: no echo at detection point left after applying the filter.');
      end
    else
      logmsg(1,'  no user-supplied coherence filters');
    end

    [idxFreeL,idxFreeR,idxPulseL,idxPulseR] = pathwayparser(pathwayList);

    % Exp.Phase contains the pulse phases in multiples of pi/2
    if ~isfield(Exp,'Phase')
      Exp.Phase = ones(1,nIntervals);  % y phase by default
    end

    % Compute all ideal pulse transfer prefactors
    pathwayprefactor = ones(1,nPathways);
    for iPulse = 1:numel(Exp.Flip)
      if idealPulse(iPulse)
        theta = Exp.Flip(iPulse)*pi/2;
        c = cos(theta/2);
        s = sin(theta/2);
        for iPathway = 1:nPathways
          switch idxPulseL(iPathway,iPulse)
            case 1, pL = c;
            case 2, pL = c;
            case 3, pL = -1i*s*(-1i)^Exp.Phase(iPulse);
            case 4, pL = -1i*s*(+1i)^Exp.Phase(iPulse);
          end
          switch idxPulseR(iPathway,iPulse)
            case 1, pR = c;
            case 2, pR = c;
            case 3, pR = +1i*s*(+1i)^Exp.Phase(iPulse);
            case 4, pR = +1i*s*(-1i)^Exp.Phase(iPulse);
          end
          pathwayprefactor(iPathway) = pL*pathwayprefactor(iPathway)*pR;
        end
      end
    end

    % Remove pathways with zero amplitude
    rmv = abs(pathwayprefactor)<1e-6;
    pathwayList(rmv,:) = [];
    idxFreeL(rmv,:) = [];
    idxFreeR(rmv,:) = [];
    idxPulseL(rmv,:) = [];
    idxPulseR(rmv,:) = [];
    pathwayprefactor(rmv) = [];

    nPathways = size(pathwayList,1);

    idxIncL = idxFreeL(:,Exp.Inc~=0);
    idxIncR = idxFreeR(:,Exp.Inc~=0);

    if EasySpinLogLevel>0
      logmsg(1,'  Pathways and prefactors:');
      Str = 'ab+-';
      for iPathway = 1:nPathways
        logmsg(1,'    %d. (%s)   %+4.3f%+4.3fi',iPathway,...
          Str(pathwayList(iPathway,:)),real(pathwayprefactor(iPathway)),imag(pathwayprefactor(iPathway)));
      end
    end

  end

  if isENDOR

    if all(Sys.lwEndor==0)
      error('Positive ENDOR line width in Sys.lwEndor needed.');
    end
    Sys.lwEndor = Sys.lwEndor(1);

    if ~isfield(Exp,'Range')
      error('Frequency range (Exp.Range) must be given for an ENDOR experiment.');
    end

    if ~isfield(Exp,'tprf')
      Exp.tprf = 20; % rf pulse length, µs
    end

  end

  % Relaxation time constants
  if ~isfield(Exp,'T1'), Exp.T1 = 0; end
  if ~isfield(Exp,'T2'), Exp.T2 = 0; end

  if predefinedExperiment
    switch Exp.ExperimentID
      case 2 % 3pESEEM
        if ~isfield(Exp,'T'), Exp.T = 0; end
      case 5 % Mims ENDOR
        if ~isfield(Exp,'T'), Exp.T = 0; end
      case 3 % 4pESEEM
        if ~isfield(Exp,'T'), Exp.T = 0; end
      case 4 % HYSCORE
        if ~isfield(Exp,'t1'), Exp.t1 = 0; end
        if ~isfield(Exp,'t2'), Exp.t2 = 0; end
        if numel(Exp.t1)>1
          error('Exp.t1 must a single positive number, in units of µs.');
        end
        if numel(Exp.t2)>1
          error('Exp.t2 must a single positive number, in units of µs.');
        end
        if Exp.t1~=Exp.t2
          fprintf('Exp.t1 and Exp.t2 are not identical, so the resulting spectrum might be asymmetric.\n');
        end
    end
  end

  OrientationSelection = isfield(Exp,'mwFreq');
  if OrientationSelection
    logmsg(1,'Microwave frequency given: orientation selection is on.');
    if ~isfield(Exp,'ExciteWidth')
      error('Orientation selection: Exp.ExciteWidth (in MHz) missing. It should be about the inverse of the first pulse length (100MHz for 10ns). If you don''t want orientation selection, set it to a very large number (1e6) or remove the microwave frequency.');
    end
  else
    logmsg(1,'no orientation selection (infinite bandwidth).');
  end

  if isfield(Exp,'ExciteWidth')
    if ~isfield(Exp,'mwFreq')
      error('Exp.ExciteWidth is given, but Exp.mwFreq is missing. Please give Exp.mwFreq.');
    end
  end

  if isfield(Exp,'HStrain')
    error('You gave Exp.HStrain, but it should be Sys.HStrain (in the system, not the experiment structure).');
  end

  %===================================================================
  % Options structure
  %===================================================================
  %

  % Nuclei: which nuclei to include in the simulation
  if isfield(Opt,'Nuclei')
    if isempty(Opt.Nuclei)
      error('Opt.Nuclei must contain the indices of nuclei to include in the simulation.');
    end
    if any(Opt.Nuclei>Sys.nNuclei) || any(Opt.Nuclei<0)
      error('Index in Opt.Nuclei out of range.');
    end
    shfNuclei = Opt.Nuclei;
  else
    % Pick nuclei to be included in the ESEEM/ENDOR computation
    shfNuclei = 1:Sys.nNuclei;
    if ~isempty(shfNuclei)
      if all(idealPulse) && isfield(Exp,'ExciteWidth')
        if Sys.fullA
          for iNuc = 1:Sys.nNuclei
            Amatrix = Sys.A((iNuc-1)*3+(1:3),:);
            idxStrongNuclei(iNuc) = max(max(abs(Amatrix)))>Exp.ExciteWidth; %#ok
          end
        else
          idxStrongNuclei = max(abs(Sys.A),[],2)>Exp.ExciteWidth;
        end
        shfNuclei(idxStrongNuclei) = [];
      end
    end
  end

  twoElectronManifolds = (Sys.nElectrons==1) && (Sys.S==1/2) && ...
    (numel(shfNuclei)==Sys.nNuclei);

  % ProductRule: determines whether product rule is used or not
  if ~isfield(Opt,'ProductRule'), Opt.ProductRule = 0; end
  if Sys.nNuclei==1, Opt.ProductRule = 0; end
  %if isENDOR, Opt.ProductRule = 1; end

  % EndorMethod: how to simulate the effect of the RF pulse in ENDOR
  if isENDOR
    if ~isfield(Opt,'EndorMethod')
      % 0 = sum-over-transitions, adjacent level population swap (wrong for >1 nucleus)
      % 1 = sum-over-transitions, bandwidth-filtered Iy pi pulse on all nuclei
      % 2 = frequency sweep, bandwidth-filtered Iy pi pulse on all nuclei
      if numel(shfNuclei)==1 || Opt.ProductRule
        Opt.EndorMethod = 0;
      else
        Opt.EndorMethod = 1;
      end
    end
    switch Opt.EndorMethod
      case 0, logmsg(1,'using population swaps of adjacent nuclear sublevels');
      case 1, logmsg(1,'using bandwidth-filtered Iy pi pulse, sum over transitions');
      case 2, logmsg(1,'using bandwidth-filtered Iy pi pulse, full RF sweep');
      otherwise, error('Unknown setting for Opt.EndorMethod. Must be 0, 1, or 2.');
    end
    if Opt.EndorMethod==0 && numel(shfNuclei)>1 && ~Opt.ProductRule
      error('Opt.EndorMethod=0 gives incorrect results for multiple nuclei.')
    end
  end

  % Expansion factor: determines size of spectral buffer
  if ~isfield(Opt,'Expand')
    if nDimensions==1, Opt.Expand = 4; else, Opt.Expand = 2; end
  end
  maxExpand = 8;
  if numel(Opt.Expand)~=1 || Opt.Expand<0 || Opt.Expand>maxExpand || rem(Opt.Expand,1)
    error('Opt.Expand must be an integer between 0 and %d.',maxExpand);
  end

  if any(realPulse) && Opt.ProductRule
    error('saffron: Cannot apply product rule and real pulses at the same time.');
  end

  if ~isfield(Opt,'OriThreshold'), Opt.OriThreshold = 0.005; end

  if ~isfield(Opt,'Window')
    Opt.Window = 'ham+';
  end

  if ~isfield(Opt,'ZeroFillFactor')
    Opt.ZeroFillFactor = 2;
  end

  processData = true;

  % undocumented options
  if ~isfield(Opt,'logplot'), Opt.logplot = 0; end
  if ~isfield(Opt,'PartialIFFT'), Opt.PartialIFFT = 1; end
  if ~isfield(Opt,'TimeDomain'), Opt.TimeDomain = 0; end


  logmsg(1,'-Hamiltonians------------------------------------------');

  %=====================================================================
  % Compute electronic Hamiltonian
  %=====================================================================

  logmsg(1,'setting up electronic Hamiltonian...');
  if twoElectronManifolds
    % not needed, S=1/2 electronic Hamiltonian can be solved analytically
  else
    coreSys = nucspinrmv(Sys,shfNuclei);
    coreSys.processed = false;
    coreSys.lw = 0;
    if isfield(coreSys,'A_')
      coreSys = rmfield(coreSys,'A');
    end
    % Operators for constructing Hamiltonian
    [H0,mux,muy,muz] = ham(coreSys);
    % Operators for computing <i|S|i>
    Sx = sop(coreSys,[1,1]); % works only for one electron spin
    Sy = sop(coreSys,[1,2]);
    Sz = sop(coreSys,[1,3]);
  end
  %=====================================================================


  %=====================================================================
  % Compute nuclear spin Hamiltonians
  %=====================================================================
  logmsg(1,'computing nuclear spin sub-Hamiltonians...');
  if ~isempty(shfNuclei)
    if Opt.ProductRule
      logmsg(1,'  separate subspace for each of the %d superhyperfine nuclei',numel(shfNuclei));
    else
      logmsg(1,'  complete nuclear state space of all %d superhyperfine nuclei',numel(shfNuclei));
      NucHams.Hnzx = 0;
      NucHams.Hnzy = 0;
      NucHams.Hnzz = 0;
      NucHams.Hhfx = 0;
      NucHams.Hhfy = 0;
      NucHams.Hhfz = 0;
      NucHams.Hnq = 0;
    end
    if isENDOR
      NucHams.Ix = 0;
      NucHams.Iy = 0;
      NucHams.Iz = 0;
    end
    I = Sys.I(shfNuclei);
    for iiNuc = numel(shfNuclei):-1:1   % only shf nuclei
      % Spin operators -------------------------------------------------
      if Opt.ProductRule
        Ix = sop(I(iiNuc),[1,1]);
        Iy = sop(I(iiNuc),[1,2]);
        Iz = sop(I(iiNuc),[1,3]);
      else
        Ix = sop(I,[iiNuc,1]);
        Iy = sop(I,[iiNuc,2]);
        Iz = sop(I,[iiNuc,3]);
      end
      iNuc = shfNuclei(iiNuc);
      % Store operators for building RF pulse
      if isENDOR
        if Opt.ProductRule
          NucHams(iiNuc).Ix = Ix;
          NucHams(iiNuc).Iy = Iy;
          NucHams(iiNuc).Iz = Iz;
        else
          NucHams.Ix = NucHams.Ix + Ix;
          NucHams.Iy = NucHams.Iy + Iy;
          NucHams.Iz = NucHams.Iz + Iz;
        end
      end

      % Nuclear Zeeman -------------------------------------------------
      pre = -Sys.gn(iNuc)*nmagn/1e3/planck/1e6; % MHz/mT
      if Opt.ProductRule
        NucHams(iiNuc).Hnzx = pre*Ix;
        NucHams(iiNuc).Hnzy = pre*Iy;
        NucHams(iiNuc).Hnzz = pre*Iz;
      else
        NucHams.Hnzx = NucHams.Hnzx + pre*Ix;
        NucHams.Hnzy = NucHams.Hnzy + pre*Iy;
        NucHams.Hnzz = NucHams.Hnzz + pre*Iz;
      end

      % Hyperfine ------------------------------------------------------
      if Sys.fullA
        A = Sys.A((iNuc-1)*3+(1:3),:);
      else
        if isfield(Sys,'AFrame')
          R_A2M = erot(Sys.AFrame(iNuc,:)).'; % A frame -> molecular frame
        else
          R_A2M = eye(3);
        end
        A = R_A2M*diag(Sys.A(iNuc,:))*R_A2M.';
        A = (A + A.')/2;
      end
      if Opt.ProductRule
        NucHams(iiNuc).Hhfx = A(1,1)*Ix + A(1,2)*Iy + A(1,3)*Iz;
        NucHams(iiNuc).Hhfy = A(2,1)*Ix + A(2,2)*Iy + A(2,3)*Iz;
        NucHams(iiNuc).Hhfz = A(3,1)*Ix + A(3,2)*Iy + A(3,3)*Iz;
      else
        NucHams.Hhfx = NucHams.Hhfx + A(1,1)*Ix + A(1,2)*Iy + A(1,3)*Iz;
        NucHams.Hhfy = NucHams.Hhfy + A(2,1)*Ix + A(2,2)*Iy + A(2,3)*Iz;
        NucHams.Hhfz = NucHams.Hhfz + A(3,1)*Ix + A(3,2)*Iy + A(3,3)*Iz;
      end

      % Nuclear quadrupole ---------------------------------------------
      Hnq_ = 0;
      if Sys.I(iNuc)>=1
        if Sys.fullQ
          Q = Sys.Q((iNuc-1)*3+(1:3),:);
        else
          if isfield(Sys,'QFrame')
            R_Q2M = erot(Sys.QFrame(iNuc,:)).'; % Q frame -> molecular frame
          else
            R_Q2M = eye(3);
          end
          Q = R_Q2M*diag(Sys.Q(iNuc,:))*R_Q2M.';
          Q = (Q + Q.')/2;
        end
        Hnq_ = ...
          Ix*(Q(1,1)*Ix + Q(1,2)*Iy + Q(1,3)*Iz) + ...
          Iy*(Q(2,1)*Ix + Q(2,2)*Iy + Q(2,3)*Iz) + ...
          Iz*(Q(3,1)*Ix + Q(3,2)*Iy + Q(3,3)*Iz);
      end
      if Opt.ProductRule
        NucHams(iiNuc).Hnq = Hnq_;
      else
        NucHams.Hnq = NucHams.Hnq + Hnq_;
      end

    end
    nSubSpaces = numel(NucHams);
    logmsg(1,'  %d nuclei, %d subspaces',numel(shfNuclei),nSubSpaces);
  else
    logmsg(1,'  no subspace factorization');
    nSubSpaces = 0;
  end
  %=====================================================================

  %=====================================================================
  % Orientation pre-selection factors
  %=====================================================================
  orientationPreSelection = OrientationSelection && twoElectronManifolds;
  if orientationPreSelection
    logmsg(1,'pre-computing orientation selecton from g tensor alone...');
    logmsg(1,'  S=1/2: using simple g tensor/HStrain model');
    logmsg(1,'  excitation width (MHz): %g',Exp.ExciteWidth);
    logmsg(1,'  HStrain (MHz): %g %g %g',Sys.HStrain(1),Sys.HStrain(2),Sys.HStrain(3));
    % g values for all orientations
    zLab = ang2vec(Orientations(:,1),Orientations(:,2));
    if Sys.fullg
      geff = Sys.g*zLab;
    else
      geff = diag(Sys.g)*zLab;
    end
    geff = sqrt(sum(geff.^2,1));
    % resonance frequencies for all orientations
    nu = geff*bmagn*Exp.Field/1e3/planck/1e6; % MHz
    % line widths for all orientations
    lw = diag(Sys.HStrain)*zLab;
    lw2 = (sum(lw.^2,1)) + Exp.ExciteWidth^2;
    gOriSelWeight = exp(-(nu-1000*Exp.mwFreq).^2./lw2);
  end

  %=====================================================================
  % Preparation
  %=====================================================================
  logmsg(1,'Preparation...');

  if isENDOR

    logmsg(1,'  ENDOR simulation');

    rf = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);
    endorspc = zeros(1,Exp.nPoints);

  else

    if Opt.TimeDomain
      logmsg(1,'  time domain simulation');
    else
      logmsg(1,'  frequency domain simulation');
      % Prepare for binning method
      ExpansionFactor = 2.^Opt.Expand;
      nPointsF = ExpansionFactor*Exp.nPoints;
      if numel(nPointsF)==2
        logmsg(1,'  %dx%d points, expand x%d -> %dx%d points',...
          Exp.nPoints(1),Exp.nPoints(2),ExpansionFactor,nPointsF(1),nPointsF(2));
      else
        logmsg(1,'  %d points, expand x%d -> %d points',Exp.nPoints,ExpansionFactor,nPointsF);
      end
      % allocate array(s)
      if nDimensions==1, siz = [1, nPointsF]; else, siz = nPointsF; end
      if Opt.ProductRule
        for iP = 1:nPathways
          for iS = 1:nSubSpaces
            pathwaybuffRe{iP,iS} = zeros(siz);
            pathwaybuffIm{iP,iS} = zeros(siz);
          end
        end
      else
        buffRe = zeros(siz);
        buffIm = zeros(siz);
      end
    end

    if nDimensions==2
      totaltd = zeros(Exp.nPoints);
    else
      totaltd = zeros(1,Exp.nPoints);
    end

  end
  %=====================================================================





  %=====================================================================
  % Orientation loop
  %=====================================================================
  logmsg(1,'Looping over %d orientations...',nOrientations);

  % Prepare offsets
  if any(realPulse)
    if ~isfield(Opt,'nOffsets')
      Opt.nOffsets = 291;
    end
    if ~isfield(Opt,'lwOffset')
      Opt.lwOffset = 100;
    end
    offsets = linspace(-1,1,Opt.nOffsets)*Opt.lwOffset*2;
    offsetWeight = exp(-(offsets/Opt.lwOffset).^2);
    offsetWeight = offsetWeight/sum(offsetWeight);
  else
    offsets = 0;
    offsetWeight = 1;
    Opt.nOffsets = 1;
  end

  if twoElectronManifolds
    g = Sys.g;
    if ~Sys.fullg
      Rg = erot(Sys.gFrame).'; % g frame -> molecular frame
      g = Rg*diag(g)*Rg.';
    end
  end

  nSkippedOrientations = 0;
  for iOri = 1:nOrientations

    % Get magnetic field orientation
    %------------------------------------------------------------------
    [~,yLab_M,zLab_M] = erot(Orientations(iOri,:),'rows');

    % Compute electronic Hamiltonian, energies and <S>
    %------------------------------------------------------------------
    if twoElectronManifolds

      if OrientationSelection
        OriSelWeight = gOriSelWeight(iOri);
      else
        OriSelWeight = 1;
      end
      if OriSelWeight<Opt.OriThreshold
        nSkippedOrientations = nSkippedOrientations + 1;
        continue
      end

      quantizationAxis_ = g.'*zLab_M;
      quantizationAxis_ = quantizationAxis_/norm(quantizationAxis_);
      Manifold(1).S = -0.5*quantizationAxis_;
      Manifold(2).S = +0.5*quantizationAxis_;

      Transitions = [1 2];
      nTransitions = 1;
      ManifoldsInvolved = [1 2];

    else

      % transition selection
      %------------------------------------------------------------
      muzL = zLab_M(1)*mux + zLab_M(2)*muy + zLab_M(3)*muz;
      H = H0 - Exp.Field*muzL;
      [eV,eE] = eig(H);
      eE = real(diag(eE));
      SyLab = yLab_M(1)*Sx + yLab_M(2)*Sy + yLab_M(3)*Sz;
      SyLab = abs(eV'*SyLab*eV);
      maxSyy = max(SyLab(:));

      if OrientationSelection
        dE = bsxfun(@minus,eE,eE.') - Exp.mwFreq*1e3; % MHz
        excitationAmplitude = exp(-(dE/Exp.ExciteWidth).^2);
        SyLab = SyLab.*excitationAmplitude;
      end

      if ~isempty(Exp.Temperature)
        % make temperature aware, include Boltzmann populations
      end

      % Remove transitions that are not excited
      SyLab(SyLab<Opt.OriThreshold*maxSyy) = 0;

      % Remove transitions that are not wanted by the user
      if ~isempty(Opt.Transitions)
        rmvTransition = ones(size(H));
        for t = 1:size(Opt.Transitions,1)
          tr = Opt.Transitions(t,:);
          rmvTransition(tr(1),tr(2)) = 0;
          rmvTransition(tr(2),tr(1)) = 0;
        end
        SyLab(rmvTransition~=0) = 0;
      end

      [v,u,OriSelWeight] = find(tril(SyLab,-1));
      Transitions = [u,v];
      nTransitions = size(Transitions,1);

      if nTransitions==0
        nSkippedOrientations = nSkippedOrientations + 1;
        continue
      end

      % computation of <S> for all manifolds involved
      %------------------------------------------------------------
      ManifoldsInvolved = zeros(1,length(Sx));
      ManifoldsInvolved(u) = 1;
      ManifoldsInvolved(v) = 1;
      ManifoldsInvolved = find(ManifoldsInvolved);
      for iM = ManifoldsInvolved
        vec = eV(:,iM);
        Manifold(iM).S = real([vec'*Sx*vec; vec'*Sy*vec; vec'*Sz*vec]);
      end
    end

    logmsg(2,'orientation %d of %d: %d transitions',iOri,nOrientations,nTransitions);

    % Compute and diagonalize nuclear Hamiltonians
    %----------------------------------------------------------------------
    for iSpace = 1:nSubSpaces
      Hnuc = Exp.Field*(zLab_M(1)*NucHams(iSpace).Hnzx + zLab_M(2)*NucHams(iSpace).Hnzy + zLab_M(3)*NucHams(iSpace).Hnzz);
      Hnuc = Hnuc + NucHams(iSpace).Hnq;
      for iM = ManifoldsInvolved
        S = Manifold(iM).S; % expectation value of spin vector, <S>
        H = Hnuc + S(1)*NucHams(iSpace).Hhfx + ...
          S(2)*NucHams(iSpace).Hhfy + ...
          S(3)*NucHams(iSpace).Hhfz;
        [VV,EE] = eig((H+H')/2);
        Manifold(iM).V{iSpace} = VV;
        Manifold(iM).E{iSpace} = real(diag(EE));
      end
    end


    % Loop over all excited EPR transitions
    %----------------------------------------------------------------------
    for iT = 1:nTransitions
      b = Transitions(iT,1); % lower manifold
      a = Transitions(iT,2); % upper manifold

      % Loop over all subspaces
      for iSpace = 1:nSubSpaces
        Ea = Manifold(a).E{iSpace};
        Eb = Manifold(b).E{iSpace};
        Ma = Manifold(a).V{iSpace};
        Mb = Manifold(b).V{iSpace};
        M = Ma'*Mb;       % <a|b> overlap matrix
        Mt = M';
        nNucStates = length(Ea);
        eyeN = eye(nNucStates);
        if any(realPulse)
          idxa = 1:nNucStates;
          idxb = idxa + nNucStates;
        end
        % Set up RF operator for ENDOR
        if isENDOR
          IyLab = yLab_M(1)*NucHams(iSpace).Ix + ...
            yLab_M(2)*NucHams(iSpace).Iy + ...
            yLab_M(3)*NucHams(iSpace).Iz;
        end

        for iOffset = 1:Opt.nOffsets

          % Pulse propagators
          for iInt = 1:nIntervals
            if realPulse(iInt)
              nu1 = (pi/2*Exp.Flip(iInt))/Exp.tp(iInt)/2/pi;
              Hpulse = [diag(Ea+offsets(iOffset)/2), +M*nu1/2i; ...
                -Mt*nu1/2i diag(Eb-offsets(iOffset)/2)];
              FullPulsePropagator = expm(-2i*pi*Exp.tp(iInt)*Hpulse);
              PulseSubPropagator{iInt,1} = FullPulsePropagator(idxa,idxa);
              PulseSubPropagator{iInt,2} = FullPulsePropagator(idxb,idxb);
              PulseSubPropagator{iInt,3} = FullPulsePropagator(idxa,idxb);
              PulseSubPropagator{iInt,4} = FullPulsePropagator(idxb,idxa);
            end
          end

          if isENDOR
            prefactor = Exp.OriWeights(iOri)*OriSelWeight(iT);
          else
            if ~Opt.ProductRule
              prefactor = Exp.OriWeights(iOri)*OriSelWeight(iT);
            else
              prefactor = 1;
            end
            if ~Opt.TimeDomain
              if Opt.ProductRule
                for iP = 1:nPathways
                  pathwaybuffRe{iP,iSpace} = zeros(siz);
                  pathwaybuffIm{iP,iSpace} = zeros(siz);
                end
              end
            end
          end
          prefactor = prefactor*offsetWeight(iOffset);

          % Compute peaks and amplitudes
          if ~predefinedExperiment

            % General method ------------------------------------------
            increments = Exp.Inc;
            BlockL = cell(0);
            BlockR = cell(0);
            for iPathway = 1:nPathways
              iBlock = 0;
              Left = prefactor*pathwayprefactor(iPathway);
              Right = 1;
              for iInt = 1:nIntervals
                % Pulse propagator
                if idealPulse(iInt)
                  switch idxPulseL(iPathway,iInt)
                    case 3, Left = M*Left;
                    case 4, Left = Mt*Left;
                  end
                  switch idxPulseR(iPathway,iInt)
                    case 3, Right = Right*Mt;
                    case 4, Right = Right*M;
                  end
                else
                  Left = PulseSubPropagator{iInt,idxPulseL(iPathway,iInt)}*Left;
                  Right = Right*PulseSubPropagator{iInt,idxPulseR(iPathway,iInt)}';
                end

                % Free evolution propagator
                t = Exp.t(iInt);
                if t>0
                  if idxFreeL(iPathway,iInt)==1, E = Ea; else, E = Eb; end
                  Left = diag(exp(-2i*pi*E*t)) * Left;
                  if idxFreeR(iPathway,iInt)==1, E = Ea; else, E = Eb; end
                  Right = Right * diag(exp(+2i*pi*E*t));
                end

                % Conclude block and start next one
                if increments(iInt)~=0
                  if iBlock==0
                    G = Left*Right;
                    if numel(G)==1
                      G = G*eyeN;
                    end
                  else
                    if numel(Left)==1
                      BlockL{iBlock} = eyeN;
                    else
                      BlockL{iBlock} = Left;
                    end
                    if numel(Right)==1
                      BlockR{iBlock} = eyeN;
                    else
                      BlockR{iBlock} = Right;
                    end
                  end
                  iBlock = iBlock + 1;
                  Left = 1;
                  Right = 1;
                end

              end % iInt

              if increments(end)~=0
                D = M;
              else
                D = Right*M*Left;
              end

              % Accumulate peaks / generate time domain
              if Opt.ProductRule
                error('Product rule for user-defined experiments not implemented.');
              end
              if Opt.TimeDomain
                totaltd = totaltd + ...
                  sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,...
                  idxIncL(iPathway,:),idxIncR(iPathway,:),...
                  Ea,Eb,G,D,BlockL{:},BlockR{:});
              else
                sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,...
                  idxIncL(iPathway,:),idxIncR(iPathway,:),...
                  Ea,Eb,G,D,BlockL{:},BlockR{:});
              end

            end % iPathway loop

          else

            switch Exp.ExperimentID

              case 5
                % Mims ENDOR ------------------------------------------------
                % coherence transfer pathway 1: +,alpha,-
                % coherence transfer pathway 2: +,beta,-

                if ~all(idealPulse)
                  error('Pre-defined Mims ENDOR with real pulses not supported.');
                end

                Q_ = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
                G_ = prefactor*Q_.*M;
                D_ = conj(Q_).*M;
                G1 = G_*Mt; D1 = D_*Mt;
                G2 = Mt*G_; D2 = Mt*D_;

                % Evolve density matrices after second pi/2 pulse
                if Exp.T~=0
                  q_ = exp(-2i*pi*Ea*Exp.T); G1 = (q_*q_').*G1;
                  q_ = exp(-2i*pi*Eb*Exp.T); G2 = (q_*q_').*G2;
                end

                % Remove nuclear coherences
                G1 = diag(diag(G1));
                G2 = diag(diag(G2));

                % Echo amplitude for off-resonant RF pulse (gives baseline)
                traceG1D1 = trace(G1*D1);
                traceG2D2 = trace(G2*D2);

                % Matrices of nuclear transition frequencies for both manifolds
                % (lower triangle is positive, upper is negative, diagonal zero)
                nu1 = abs(bsxfun(@minus,Ea,Ea.'));
                nu2 = abs(bsxfun(@minus,Eb,Eb.'));

                % Transform RF pulse Hamiltonians from Zeeman to eigenbasis
                Iy1 = Ma'*IyLab*Ma;
                Iy2 = Mb'*IyLab*Mb;

                fwhm = 1/Exp.tprf; % excitation bandwidth, MHz
                Gamma = fwhm/(2*sqrt(log(2)));
                theta = pi; % rf pulse flip angle

                Gam = Sys.lwEndor/sqrt(2*log(2));  % line broadening
                pre = sqrt(2/pi)/Gam;

                switch Opt.EndorMethod
                  case 2
                    % Sweep method
                    %---------------------------------------------------------
                    % Sweep rf, apply bandwidth-limited rf pulse at each
                    % rf frequency that is within bandwidth of a nuclear
                    % transition frequency.
                    BWthreshold = 0.01;
                    for irf = 1:numel(rf) % do full rf sweep
                      % calculate bandpass filter matrices (Gaussian profile)
                      BW1 = exp(-((nu1-rf(irf))/Gamma).^2);
                      BW2 = exp(-((nu2-rf(irf))/Gamma).^2);
                      if any(BW1(:)>BWthreshold)
                        Prfa = expm(-1i*theta*(Iy1.*BW1));
                        G1_ = Prfa*G1*Prfa'; % apply rf pulse
                        G1_ = diag(diag(G1_)); % remove nuclear coherences
                        endorspc(irf) = endorspc(irf) + traceG1D1 - trace(G1_*D1);
                      end
                      if any(BW2(:)>BWthreshold)
                        Prfb = expm(-1i*theta*(Iy2.*BW2));
                        G2_ = Prfb*G2*Prfb'; % apply rf pulse
                        G2_ = diag(diag(G2_)); % remove nuclear coherences
                        endorspc(irf) = endorspc(irf) + traceG2D2 - trace(G2_*D2);
                      end
                    end

                  case 0
                    % Sum-over-transitions method, adjacent-level population swap
                    %------------------------------------------------------------
                    % Loop over all pairs of adjacent nuclear levels and swap
                    % populations. This is correct only for a single nucleus at a
                    % time. (only available method prior to 5.0.21)

                    % loop only over adjacent nuclear sublevel pairs
                    for j = 1:nNucStates-1
                      for i = j+1
                        % RF pulse approximation: only swap diagonal elements ii and jj
                        G1_ = G1; q = G1_(i,i); G1_(i,i) = G1_(j,j); G1_(j,j) = q;
                        G2_ = G2; q = G2_(i,i); G2_(i,i) = G2_(j,j); G2_(j,j) = q;
                        ampl = [traceG1D1-trace(G1_*D1), traceG2D2-trace(G2_*D2)];
                        freq = [Ea(i)-Ea(j), Eb(i)-Eb(j)];
                        %idx = fix(Exp.nPoints*(freq-rf(1))/(rf(end)-rf(1)))+1;
                        %endorspc(idx) = endorspc(idx) + ampl;
                        %endorspc = endorspc + ...
                        %  lisum1i(Template.y,Template.x0,Template.lw,freq,ampl,Sys.lwEndor*[1 1],rf);
                        endorspc = endorspc + ampl(1)*pre*exp(-2*((rf-freq(1))/Gam).^2);
                        endorspc = endorspc + ampl(2)*pre*exp(-2*((rf-freq(2))/Gam).^2);
                      end
                    end

                  case 1
                    % Sum-over-transitions method, bandwidth-limited Iy pulse
                    %---------------------------------------------------------
                    % Loop over all nuclear transition and apply
                    % bandwidth-limited rf pulse operator. The excitation
                    % bandwidth is Gaussian.
                    [i,j] = find(triu(ones(nNucStates),1));
                    for iNucTrans = 1:numel(i)
                      % Calculate propagators for narrow-band Gaussian
                      % pulses at the nuclear transition frequency i->j
                      % for both manifolds
                      freq1 = nu1(i(iNucTrans),j(iNucTrans));
                      freq2 = nu2(i(iNucTrans),j(iNucTrans));
                      BW1 = exp(-((nu1-freq1)/Gamma).^2);
                      BW2 = exp(-((nu2-freq2)/Gamma).^2);
                      Prfa = expm(-1i*theta*(Iy1.*BW1));
                      Prfb = expm(-1i*theta*(Iy2.*BW2));
                      % Propagate densities and remove all coherences
                      G1_ = diag(diag(Prfa*G1*Prfa'));
                      G2_ = diag(diag(Prfb*G2*Prfb'));
                      % Calculate ENDOR peak amplitudes
                      sig1 = traceG1D1 - trace(G1_*D1);
                      sig2 = traceG2D2 - trace(G2_*D2);
                      %endorspc = endorspc + ...
                      %  lisum1i(Template.y,Template.x0,Template.lw,freq,ampl,Sys.lwEndor*[1 1],rf);
                      % Add peaks to spectrum
                      endorspc = endorspc + sig1*pre*exp(-((rf-freq1)/Gam).^2);
                      endorspc = endorspc + sig2*pre*exp(-((rf-freq2)/Gam).^2);
                    end
                end

              case 4
                % HYSCORE ---------------------------------------------------
                % coherence transfer pathway 1: +,alpha,beta,-
                % coherence transfer pathway 2: +,beta,alpha,-

                if ~all(idealPulse)
                  error('Pre-defined HYSCORE with real pulses not supported.');
                end

                Q_ = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
                G_ = prefactor*Q_.*M;
                D_ = conj(Q_).*M;
                G1 = G_*Mt; G2 = Mt*G_;
                D1 = Mt*D_; D2 = D_*Mt;
                Tl1 = Mt; Tr1 = M;
                Tl2 = M; Tr2 = Mt;

                if Exp.t1~=0
                  q_ = exp(-2i*pi*Exp.t1*Ea); G1 = (q_*q_').*G1;
                  q_ = exp(-2i*pi*Exp.t1*Eb); G2 = (q_*q_').*G2;
                end
                if Exp.t2~=0
                  q_ = exp(+2i*pi*Exp.t2*Eb); D1 = (q_*q_').*D1;
                  q_ = exp(+2i*pi*Exp.t2*Ea); D2 = (q_*q_').*D2;
                end

                if Opt.TimeDomain
                  if Opt.ProductRule
                    pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    pathwaytd{2,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  else
                    totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  end
                else
                  if Opt.ProductRule
                    sf_peaks(IncSchemeID,pathwaybuffRe{1,iSpace},pathwaybuffIm{1,iSpace},Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    sf_peaks(IncSchemeID,pathwaybuffRe{2,iSpace},pathwaybuffIm{2,iSpace},Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  else
                    sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  end
                end

              case 3
                % 4pESEEM ---------------------------------------------------
                % coherence transfer pathway 1: +,alpha,beta,-
                % coherence transfer pathway 2: +,beta,alpha,-

                if ~all(idealPulse)
                  error('Pre-defined 4p-ESEEM with real pulses not supported.');
                end

                Q_ = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
                G_ = prefactor*Q_.*M;
                D_ = conj(Q_).*M;
                G1 = G_*Mt; D1 = Mt*D_;
                G2 = Mt*G_; D2 = D_*Mt;
                Tl1 = Mt; Tr1 = M;
                Tl2 = M; Tr2 = Mt;

                if Exp.T~=0
                  q_ = exp(-2i*pi*Exp.T*Ea); G1 = (q_*q_').*G1;
                  q_ = exp(-2i*pi*Exp.T*Eb); G2 = (q_*q_').*G2;
                  q_ = exp(+2i*pi*Exp.T*Eb); D1 = (q_*q_').*D1;
                  q_ = exp(+2i*pi*Exp.T*Ea); D2 = (q_*q_').*D2;
                end

                if Opt.TimeDomain
                  if Opt.ProductRule
                    pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    pathwaytd{2,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  else
                    totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  end
                else
                  if Opt.ProductRule
                    sf_peaks(IncSchemeID,pathwaybuffRe{1,iSpace},pathwaybuffIm{1,iSpace},Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    sf_peaks(IncSchemeID,pathwaybuffRe{2,iSpace},pathwaybuffIm{2,iSpace},Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  else
                    sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,[1 2],[1 2],Ea,Eb,G1,D1,Tl1,Tr1);
                    sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,[2 1],[2 1],Ea,Eb,G2,D2,Tl2,Tr2);
                  end
                end

              case 2
                % 3pESEEM ------------------------------------------------------------
                % coherence transfer pathway 1: +,alpha,-
                % coherence transfer pathway 2: +,beta,-

                if ~all(idealPulse)
                  error('Pre-defined 3pESEEM not supported for real pulses.');
                end

                Q_ = exp(-2i*pi*Ea*Exp.tau)*exp(-2i*pi*Eb*Exp.tau)';
                G_ = prefactor*Q_.*M;
                D_ = conj(Q_).*M;
                G1 = G_*Mt; D1 = D_*Mt;
                G2 = Mt*G_; D2 = Mt*D_;

                if Exp.T~=0
                  q_ = exp(-2i*pi*Exp.T*Ea); G1 = (q_*q_').*G1;
                  q_ = exp(-2i*pi*Exp.T*Eb); G2 = (q_*q_').*G2;
                end

                if Opt.TimeDomain
                  if Opt.ProductRule
                    pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,1,1,Ea,Eb,G1,D1);
                    pathwaytd{2,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,2,2,Ea,Eb,G2,D2);
                  else
                    totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,1,1,Ea,Eb,G1,D1);
                    totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,2,2,Ea,Eb,G2,D2);
                  end
                else
                  if Opt.ProductRule
                    sf_peaks(IncSchemeID,pathwaybuffRe{1,iSpace},pathwaybuffIm{1,iSpace},Exp.dt,1,1,Ea,Eb,G1,D1);
                    sf_peaks(IncSchemeID,pathwaybuffRe{2,iSpace},pathwaybuffIm{2,iSpace},Exp.dt,2,2,Ea,Eb,G2,D2);
                  else
                    sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,1,1,Ea,Eb,G1,D1);
                    sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,2,2,Ea,Eb,G2,D2);
                  end
                end

              case 1
                % 2pESEEM ---------------------------------------------------------
                % coherence transfer pathway: +-

                if ~all(idealPulse)
                  error('Pre-defined 2pESEEM not supported for real pulses.');
                end
                G = prefactor*M;
                D = M;
                T1left = Mt;
                T1right = Mt;
                if Exp.tau>0
                  Q_ = exp(-2i*pi*Exp.tau*Ea)*exp(-2i*pi*Exp.tau*Eb)';
                  G = Q_.*G;
                  D = conj(Q_).*D;
                end

                if Opt.TimeDomain
                  if Opt.ProductRule
                    pathwaytd{1,iSpace} = sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                  else
                    totaltd = totaltd + sf_evolve(IncSchemeID,Exp.nPoints,Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                  end
                else
                  if Opt.ProductRule
                    sf_peaks(IncSchemeID,pathwaybuffRe{1,iSpace},pathwaybuffIm{1,iSpace},Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                  else
                    sf_peaks(IncSchemeID,buffRe,buffIm,Exp.dt,[1 2],[2 1],Ea,Eb,G,D,T1left,T1right);
                  end
                end

            end % Exp.ExperimentID switchyard

          end % if predefinedExperiment
        end % offset loop
      end % nuclear subspace loop


      % Apply product rule
      %--------------------------------------------------------------
      if Opt.ProductRule
        if ~isENDOR

          for iPathway = 1:nPathways
            td_ = 1;
            for iSpace = 1:nSubSpaces
              if Opt.TimeDomain
                thistd = pathwaytd{iPathway,iSpace};
              else
                pathwaybuff = complex(pathwaybuffRe{iPathway,iSpace},pathwaybuffIm{iPathway,iSpace});
                if nDimensions==1
                  thistd = ifft(pathwaybuff)*nPointsF;
                  thistd = thistd(1:Exp.nPoints);
                else
                  thistd = ifft2dpartial(pathwaybuff,Exp.nPoints,Opt.PartialIFFT)*prod(nPointsF);
                  %thistd = thistd;
                end
              end
              td_ = td_.*thistd;
            end
            totaltd = totaltd + Exp.OriWeights(iOri)*OriSelWeight(iT)*td_;
          end

        else
          % ENDOR w. product rule: no action here
        end % if ~isENDOR

      else % no product rule
        % TD: nothing
        % FD: IFFT is done after orientation loop
      end % if Opt.ProductRule

    end % electronic transition loop

  end % orientation loop

  logmsg(1,'end of orientation/transition loop');
  logmsg(1,'%d of %d orientations skipped',nSkippedOrientations,nOrientations);
  %=================================================================


  %=================================================================
  % Postprocessing
  %=================================================================
  if isENDOR

    endorspc = real(endorspc);

    % Normalize modulation signal
    % No need to normalize out the experiment prefactors (due to
    % pulse transfer amplitudes), since they were not included above.
    EqDensityTrace = prod(2*Sys.I+1);
    endorspc = endorspc/EqDensityTrace;
    endorspc = endorspc/nPathways;

    if Opt.EndorMethod==2
      % convolve only for sweep method (not needed for sum-over-transitions)
      endorspc = convspec(endorspc,rf(2)-rf(1),Sys.lwEndor);
    end

  else

    if Opt.ProductRule
      td = totaltd;
    else
      if Opt.TimeDomain
        td = totaltd;
      else
        logmsg(1,'Postprocessing...');
        buff = complex(buffRe,buffIm);
        if nDimensions==1
          td = ifft(buff)*numel(buff);
          td = td(1:Exp.nPoints);
        else
          td = ifft2dpartial(buff,Exp.nPoints,Opt.PartialIFFT)*numel(buff);
        end
      end
    end

    % Normalize modulation signal
    % No need to normalize out the experiment prefactors (due to
    % pulse transfer amplitudes), since they were not included above.
    EqDensityTrace = prod(2*Sys.I+1);
    td = td/EqDensityTrace;
    td = td/nPathways;

    if ~predefinedExperiment
      td = td/max(abs(pathwayprefactor));
    end

  end

  %=================================================================
  % Time axes, relaxation
  %=================================================================
  if ~isENDOR
    switch Exp.ExperimentID
      case -1
        switch nDimensions
          case 1
            t1 = (0:Exp.nPoints-1)*Exp.dt;
          case 2
            t1 = (0:Exp.nPoints(1)-1)*Exp.dt(1);
            t2 = (0:Exp.nPoints(2)-1)*Exp.dt(2);
        end
      case 1 % 2p-ESEEM
        t1 = (0:Exp.nPoints-1)*Exp.dt + Exp.tau(1);
      case 2 % 3p-ESEEM
        t1 = (0:Exp.nPoints-1)*Exp.dt + Exp.tau(1) + Exp.T;
      case 3 % 4p-ESEEM
        t1 = (0:Exp.nPoints-1)*Exp.dt + Exp.T;
      case 4 % HYSCORE
        t1 = (0:Exp.nPoints(1)-1)*Exp.dt(1) + Exp.t1;
        t2 = (0:Exp.nPoints(2)-1)*Exp.dt(2) + Exp.t2;
      case 5 % Mims ENDOR
        t1 = (0:Exp.nPoints-1)*Exp.dt;
    end

    decayAdded = false;
    if any(~isinf(Sys.T1T2))
      logmsg(1,'Adding relaxation decays...');
      T1 = Sys.T1T2(1);
      T2 = Sys.T1T2(2);
      tdecay = [];
      switch Exp.ExperimentID
        case 1 % two-pulse ESEEM
          if ~isinf(T2)
            tdecay = exp(-2*t1/T2);
            td = td.*tdecay;
          end
        case 2 % three-pulse ESEEM
          if ~isinf(T1)
            tdecay = exp(-2*Exp.tau/T2)*exp(-t1/T1);
            td = td.*tdecay;
          end
        case 4 % HYSCORE
          if ~isinf(T1)
            tdecay = exp(-t1/T1)'*exp(-t2/T1);
            tdecay = exp(-2*Exp.tau/T2)*tdecay;
            td = td.*tdecay;
          end
      end
      decayAdded = ~isempty(tdecay);
    end

  end


  %===============================================================
  % TD data processing
  %===============================================================
  logmsg(1,'-final-------------------------------------------------');
  logmsg(1,'Data processing...');
  if ~isENDOR
    switch nDimensions
      case 1
        if processData

          % Decay correction
          if decayAdded
            [~,~,tdfit] = exponfit(t1,td,1);
            tdx = td - tdfit;
          else
            tdx = td;
          end

          % Baseline correction
          tdx = tdx - mean(tdx);
          %plot(t1,tdx); pause;
          %tdx(1) = tdx(1)/2;

          % Apodization
          win = apowin(Opt.Window,numel(tdx)).';
          tdx = tdx.*win;

          % Fourier transformation
          fd = fft(tdx,Opt.ZeroFillFactor*numel(tdx));
          fd = fftshift(fd);
          f1 = fdaxis(Exp.dt,length(fd));
        end

      case 2
        if processData
          if decayAdded
            tdx = basecorr(td,[1 2],[2 2]);
          else
            tdx = basecorr(td,[1 2],[0 0]);
          end

          w1 = apowin(Opt.Window,Exp.nPoints(1));
          w2 = apowin(Opt.Window,Exp.nPoints(2));

          fd = fftshift(fftn(tdx.*(w1*w2.'),Opt.ZeroFillFactor*Exp.nPoints));
          f1 = fdaxis(Exp.dt(1),size(fd,1));
          f2 = fdaxis(Exp.dt(2),size(fd,2));
        end

    end

    if max(abs(fd))<1e-300
      fd = fd*0;
    end

    % Collect output structure
    if processData
      if nDimensions==2, info.f1 = f1; info.f2 = f2; else, info.f = f1; end
      info.fd = fd;
    else
      info = struct;
    end

  else
    % Collect output structure
    f1 = rf;
    fd = endorspc;
    if max(abs(fd))<1e-300
      fd = fd*0;
    end
    info = struct;
  end

  % Assign output
  %---------------------------------------------------------------
  if isENDOR
    x_out = f1;
    y_out = fd;
  else
    if nDimensions==2, x_out = {t1,t2}; else, x_out = t1; end
    y_out = td;
  end

  % Output
  if isENDOR
    info.fd = fd;
    info.td = [];
  else
    info.td = td;
    info.fd = fd;
  end

  switch nargout
    case 1, varargout = {y_out};
    case 2, varargout = {x_out,y_out};
    case 3, varargout = {x_out,y_out,info};
  end

  %=============================================================================
else  % if fastSimulationMode

  logmsg(1,'-processing Sys structure------------------------------');

  [Sys, Sigma, DetOps, Events, Relaxation] = s_propagationsetup(Sys,Events,Opt);

  % Prepare for the frame shift -----------------------------------------
  nElectrons = length(Sys.S);

  if Opt.FrameShift ~= 0
    logmsg(1,'  adapting Sys.g to the simulation frame');
  end
  gshift = (Opt.FrameShift*1e9)*planck/bmagn/(Exp.Field(end)*1e-3);

  issize = @(A,siz) isequal(size(A),siz);
  fullg = issize(Sys.g,[3*nElectrons 3]);
  if fullg
    gshiftmat = repmat(gshift*eye(3),[nElectrons,1]);
    Sys.g = Sys.g - gshiftmat;
  else
    Sys.g = Sys.g - gshift;
  end
  %----------------------------------------------------------------------

  rawSignals = cell(1,nOrientations);
  timeAxis = cell(1,nOrientations);

  logmsg(1,'-starting orientation loop-----------------------------');
  if Opt.Relaxation
    logmsg(1,'  relaxation is active during simulation');
  else
    logmsg(1,'  no relaxation during propagation');
  end
  logmsg(1,'  propagating %d orientations',nOrientations);
  Field = Exp.Field;

  parfor iOrientation = 1 : nOrientations

    Sys_L = rotatesystem(Sys,Orientations(iOrientation,:)); % mol->lab frame

    Ham = ham(Sys_L,Field*[0 0 1]);  % lab frame

    Relaxation_ = Relaxation;
    if ~isempty(Relaxation_)
      logmsg(1,'  adapting relaxation superoperator to system frame');
      [U,~] = eig(Ham);
      R = kron(transpose(U),U');
      Relaxation_.Gamma = R'*Relaxation_.Gamma*R;
    end

    [timeAxis{iOrientation}, rawSignals{iOrientation}] = s_thyme(Sigma, Ham, DetOps, Events, Relaxation_, Vary);
  end

  timeAxis = timeAxis{1};

  % Accumulate powder signal
  Signal = 0;
  for iOrientation = 1 : nOrientations
    Signal = Signal + rawSignals{iOrientation}*Exp.OriWeights(iOrientation);
  end

  if isfield(Exp,'DetPhase')
    logmsg(1,'  applying detection phase: %d*pi',Exp.DetPhase/pi);
    phase = exp(-1i*Exp.DetPhase);
  else
    logmsg(1,'  applying default detection phase: 0');
    phase = exp(-1i*0);
  end

  Signal = Signal*phase;

  if ~Opt.SinglePointDetection && isfield(Exp,'DetIntegrate') && Exp.DetIntegrate
    logmsg(1,'-integrating echos-------------------------------------');

    sizeSig = size(Signal);

    Signal = reshape(Signal,[nDataPoints, sizeSig(end)]);

    Signal = sum(Signal,2);

    Signal = reshape(Signal, [sizeSig(1:end-2) 1]);

    Opt.SinglePointDetection = true;

    logmsg(1,'  integrated %d transients',numel(Signal));
  end

  % Signal processing
  if ~Opt.SinglePointDetection
    logmsg(1,'-processing transients---------------------------------');

    if ~isfield(Exp,'DetFrequency') && isfield(Exp,'mwFreq')
      logmsg(1,'  using Exp.mwFreq as detection frequency');
      Exp.DetFrequency = Exp.mwFreq;
    end

    % Adapt FreqTranslation if needed
    FreqTranslation = -(Exp.DetFrequency - Events{1}.FrameShift);

    % Downconversion/processing of signal
    Signal = signalprocessing(timeAxis,Signal,FreqTranslation);

    %       Signal = Signal/max(abs(Signal(:)));
    if ndims(Signal) == 2 && size(Signal,2) ==1 %#ok<ISMAT>
      Signal = permute(Signal,[2 1]);
    end

    dt = Events{end}.TimeStep;
    timeAxis = Exp.DetWindow(1):dt:Exp.DetWindow(2);
  end

  logmsg(1,'-setting up axes for output----------------------------');

  if ~isfield(Exp,'nPoints')
    logmsg(1,'  no indirect dimensions');
    nIndirectDimensions = 0;
    x = [];  % one datapoint
  else
    nIndirectDimensions = length(Exp.nPoints);
    logmsg(1,'  setting up %d axes for indirect dimensions',nIndirectDimensions);
    x = cell(1,nIndirectDimensions);
    predefinedExperiment = isfield(Exp,'Sequence') && ~isempty(Exp.Sequence) && ischar(Exp.Sequence);
    if predefinedExperiment
      ExperimentNames = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR'};
      Exp.ExperimentID = find(strcmp(Exp.Sequence,ExperimentNames));
      switch Exp.ExperimentID
        case 1  % 2pESEEM
          if isfield(Exp,'tau')
            taustart = Exp.tau;
          else
            taustart = 0;
          end
          x{1} = taustart + Exp.dt*(0:Exp.nPoints-1);
        case 2  % 3pESEEM
          if isfield(Exp,'T')
            Tstart = Exp.T;
          else
            Tstart = 0;
          end
          if isfield(Exp,'tau')
            taustart = Exp.tau;
          else
            taustart = 0;
          end
          x{1} = taustart + Tstart + Exp.dt*(0:Exp.nPoints-1);
        case 3  % 4pESEEM
          if isfield(Exp,'T')
            Tstart = Exp.T;
          else
            Tstart = 0;
          end
          x{1} = Tstart + Exp.dt*(0:Exp.nPoints-1);
        case 4  % HYSCORE
          for iDim = 1:2
            if isfield(Exp,strcat('t',num2str(iDim)))
              tstart = Exp.(strcat('t',num2str(iDim)));
            else
              tstart = 0;
            end
            x{iDim} = tstart + Exp.dt*(0:Exp.nPoints(iDim)-1);
          end
      end
    else
      for iDim = 1:nIndirectDimensions
        Dimension_ = ['Dim' num2str(iDim)];
        nParameters = size(Exp.(Dimension_),1);
        axes_ = zeros(nParameters,Exp.nPoints(iDim));
        for iParameter = 1: nParameters
          if length(Exp.(Dimension_){iParameter,2}) == 1
            axes_(iParameter,:) = Exp.(Dimension_){iParameter,2}*(0:Exp.nPoints(iDim)-1);
          else
            axes_(iParameter,:) = 1:Exp.nPoints(iDim);
          end
        end
        x{iDim} = axes_;
      end
    end
  end

  % add the transient time axis to the output
  if ~Opt.SinglePointDetection
    logmsg(1,'  adding axis of direct dimension');
    x{nIndirectDimensions+1} = timeAxis;
  end

  if iscell(x) && numel(x)==1
    x = x{1};
  end

  info = struct;
  info.SinglePointDetection = Opt.SinglePointDetection;

  switch nargout
    case 0, s_plotting(timeAxis,Signal,Exp,Opt);
    case 1, varargout = {Signal};
    case 2, varargout = {x,Signal};
    case 3, varargout = {x,Signal,info};
  end

end

%===============================================================
% Report performance
%===============================================================
endTime = datetime;
elapsedtime = endTime-startTime;
logmsg(1,'saffron took %s',elapsedtime);

logmsg(1,'=end=saffron======%s=================\n',datetime);

clear global EasySpinLogLevel

end
%===============================================================================
%===============================================================================
%===============================================================================


%=======================================================================
% Performs a 2D Inverse Fourier Transformation of a 2D array spc and
% returns the (1:N,1:N) portion of the result. Efficient when N is
% much smaller than the array dimension.
%=======================================================================
% is faster for (128*4, 128*8, 128*16; 256*4, 256*8; 512*4, 512*8)
function td = ifft2dpartial(spc,N,usePartialIFFT)
if usePartialIFFT
  % IFFT along one dim, pick 1:N, IFFT along the other, pick 1:N
  td = ifft(spc,[],1);
  td = ifft(td(1:N(1),:),[],2);
  td = td(:,1:N(2));
else
  % IFFT2 the entire array and then pick small block (1:N,1:N)
  td = ifft2(spc);
  td = td(1:N(1),1:N(2));
end
end

%==========================================================================
%==========================================================================
%==========================================================================
function Signal = sf_evolve(IncSchemeID,nPoints,dt,incL,incR,Ea,Eb,G,D,T1l,T1r)

if nargin-9~=2*(numel(incL)-1)
  error('Time-domain evolution kernel: Wrong number of transfer matrices.');
end

NN = numel(G);
Density = G;
Detector = reshape(D.',1,NN);
E = {Ea,Eb};

% Pre-allocate signal array
if numel(nPoints)==1
  Signal = zeros(1,nPoints);
else
  Signal = zeros(nPoints);
end

switch IncSchemeID

  case 1 % IncScheme [1], three-pulse ESEEM, five-pulse ESEEM
    FinalDensity = Density(:);
    UUt = exp(-2i*pi*E{incL(1)}*dt)*exp(-2i*pi*E{incR(1)}*dt)';
    UUt = UUt(:);
    for k1 = 1:nPoints
      Signal(k1) = Detector*FinalDensity;
      FinalDensity = UUt.*FinalDensity;
    end

  case 2 % IncScheme [1 1], two-pulse ESEEM, 1D-CP, refocused two-pulse ESEEM
    UUleft  = exp(-2i*pi*E{incL(2)}*dt)*exp(-2i*pi*E{incL(1)}*dt).';
    UUright = exp(+2i*pi*E{incR(1)}*dt)*exp(+2i*pi*E{incR(2)}*dt).';
    for k = 1:nPoints
      % compute density right before detection
      FinalDensity = T1l*Density*T1r;
      % compute trace(Detector*FinalDensity)
      Signal(k) = Detector*FinalDensity(:);
      % Evolve transfer matrices forward and backwards at the same time.
      T1l = UUleft.*T1l; % equivalent to U*T1Left*U'
      T1r = UUright.*T1r; % equivalent to U*T1Right*U'
    end

  case 3 % IncScheme [1 -1]
    error('Time-domain evolution for incrementation scheme [1 -1] not supported.');

  case 11 % IncScheme [1 2], HYSCORE
    UUt1 = exp(-2i*pi*E{incL(1)}*dt(1))*exp(-2i*pi*E{incR(1)}*dt(1))';
    UUt2 = exp(-2i*pi*E{incL(2)}*dt(2))*exp(-2i*pi*E{incR(2)}*dt(2))';
    UUt2 = UUt2(:);
    for k1 = 1:nPoints(1)
      FinalDensity = reshape(T1l*Density*T1r,NN,1);
      for k2 = 1:nPoints(2)
        Signal(k1,k2) = Detector*FinalDensity;
        FinalDensity = UUt2.*FinalDensity;
      end
      Density = UUt1.*Density;
    end

  case 12 % IncScheme [1 2 2]
    error('Time-domain evolution for this experiment not supported.');

  case 13 % IncScheme [1 2 1] 2D-3p
    error('Time-domain evolution for this experiment not supported.');

  case 14 % IncScheme [1 1 2], CF-NF
    error('Time-domain evolution for this experiment not supported.');

  case 15 % IncScheme [1 2 2 1], 2D-CP
    error('Time-domain evolution for this experiment not supported.');

  case 16 % IncScheme [1 2 -2 1]
    error('Time-domain evolution for this experiment not supported.');

  case 17 % IncScheme [1 1 2 2], 2D refocused primary ESEEM
    error('Time-domain evolution for this experiment not supported.');

  otherwise

    error('Time-domain evolution for this experiment not supported.');

end

end


%-------------------------------------------------------------------------------
% Coherence transfer pathway parser
%-------------------------------------------------------------------------------
function [FreeL,FreeR,PulsL,PulsR] = pathwayparser(ctp)

% ctp: a 2D char array with one CTP per row
%   two-pulse ESEEM: '+-'
%   three-pulse ESEEM: ['+a-';'+b-']
%   HYSCORE: ['+ab-';'+ba'] or ['+ab-';'+ba-';'+aa-';'+bb-']
%   etc.

% General superpropagator
%      Sa  Sb  S+  S-
% Sa   aa  ++  a+  +a
% Sb   --  bb  -b  b-
% S+   a-  +b  ab  +-
% S-   -a  b+  -+  ba
%
% E.g. -b in 2nd row and 3rd column describes transfer from
% S+ to Sb. -b is an abbreviation for U- * ... * Ub'.

Display = 0;
idealPulse = nargout<=2;

% Pre- and post-multiplications
PropL = [1 3 1 3; 4 2 4 2; 1 3 1 3; 4 2 4 2];
PropR = [1 3 3 1; 4 2 2 4; 4 2 2 4; 1 3 3 1];

% Convert ctp string with ab+- to number array with 1234
% General code:
% alpha = a = 1, beta  = b = 2, plus  = + = 3, minus = - = 4
if ischar(ctp)
  code0('ab+-') = [1 2 3 4];
  ctp = code0(ctp);
end

% Get indices for pre- and post-multiplication for pulses
if ~idealPulse
  before = [ones(size(ctp,1),1)*2 ctp(:,1:end-1)];
  idx = ctp + (before-1)*4;
  PulsL = PropL(idx);
  PulsR = PropR(idx);
end

% Get indices for pre- and post-multiplications for free evolutions
idx = 5*ctp - 4;
FreeL = PropL(idx);
FreeR = PropR(idx);

% Display formula
if Display
  decode = 'ab+-';
  ecode = 'AB';
  for i = 1:size(ctp,1)
    formula = '';
    for k = 1:size(ctp,2)
      if ~Opt.idealPulse
        formula = ['P' decode(PulsL(i,k)) ' ' formula ' P' decode(PulsR(i,k)) ''''];
      end
      formula = [ecode(FreeL(i,k)) ' ' formula ' ' ecode(FreeR(i,k)) ''''];
    end
    fprintf('%s :   trace(%s)\n',decode(ctp(i,:)),formula);
  end
end

end


% Returns spin system Sys transformed from the molecular frame to the lab
% frame, with Euler angles angles_M2L = [phi theta chi] that describe the
% transformation from molecular frame to lab frame.
function Sys_L = rotatesystem(Sys_M,angles_M2L)

Sys_L = Sys_M;

R_M2L = erot(angles_M2L);  % molecular frame -> lab frame
R_L2M = R_M2L.'; % lab frame -> molecular frame

tensors = {'g', 'A', 'Q', 'D', 'ee', 'nn'};
tensors2rotate = tensors(isfield(Sys_M,tensors));

for t = 1:numel(tensors2rotate)
  tFrame_ = [tensors2rotate{t} 'Frame'];
  ang_M2T = Sys_M.(tFrame_);
  for iRow = 1:size(ang_M2T,1)
    for iCol = 1:size(ang_M2T,2)/3
      colidx = 3*(iCol-1)+1:3*iCol;
      R_M2T = erot(ang_M2T(iRow,colidx)); % molecular frame -> tensor eigenframe
      R_L2T = R_M2T*R_L2M;  % lab frame -> tensor eigenframe
      ang_L2T(iRow,colidx) = eulang(R_L2T);
    end
  end
  Sys_L.(tFrame_) = ang_L2T;
end

end

%===============================================================================
% Plotting function
%===============================================================================
function saffron_plot(x,info,isENDOR,Sys,Exp,Opt)

logmsg(1,'Plotting...');

twoDim = iscell(x);
if twoDim
  x1 = x{1};
  x2 = x{2};
else
  x1 = x;
end

clf

if isENDOR

  plot(x1,info.fd);
  xlabel('frequency (MHz)');
  ylabel('intensity (arb.u.)');
  if isfield(Exp,'tau') && (Exp.tau<1)
    h=title(sprintf('Mims ENDOR, %g mT, tau = %g ns',Exp.Field,1000*Exp.tau));
  else
    h=title(sprintf('Mims ENDOR, %g mT, tau = %g µs',Exp.Field,Exp.tau));
  end
  set(h,'FontWeight','b');

else

  plotFreqDomain = (isfield(info,'f') || isfield(info,'f1')) && isfield(info,'fd');

  if ~twoDim

    % Time 
    if plotFreqDomain
      subplot(2,1,1);
    end
    predefinedExperiment = isfield(Exp,'ExperimentID') && Exp.ExperimentID>0;
    ExperimentNames = {'2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR'};
    plotQuadratureSignal = ~predefinedExperiment && ~isreal(info.td);
    if plotQuadratureSignal
      plot(x1,real(info.td),x1,imag(info.td));
      legend('I','Q');
      legend boxoff
    else
      plot(x1,real(info.td));
    end
    axis tight
    xl = xlim;
    yl = ylim;
    ylim(yl+[-1 1]*diff(yl)*0.1);

    if predefinedExperiment
      xlim([0 xl(2)]);
      xlb = {'\tau (µs)','\tau+T (µs)','T (µs)','...','frequency (MHz)'};
      xlabel(xlb{Exp.ExperimentID});
      ylabel('echo amplitude');
      title([ExperimentNames{Exp.ExperimentID},', TD signal']);
    else
      xlim([x1(1) x1(end)])
      xlabel('t (µs)');
      ylabel('echo amplitude (arb.u.)');
      title('User-defined experiment, TD signal');
    end
    set(gca,'Layer','top');

    % Frequency domain
    if plotFreqDomain
      subplot(2,1,2);
      idx = find(info.f==0):length(info.f);
      xf = info.f(idx);
      if plotQuadratureSignal
        h = plot(xf,abs(info.fd(idx)),xf,real(info.fd(idx)),xf,imag(info.fd(idx)));
        legend('abs','in-phase','quadrature');
        legend boxoff
      else
        h = plot(xf,abs(info.fd(idx)),xf,real(info.fd(idx)));
        legend('abs','in-phase');
        legend boxoff
      end
      axis tight
      xlim([0 max(info.f)]);
      xlabel('\nu (MHz)');
      ylabel('intensity (arb.u.)');
      title('Spectrum');
      if isfield(Sys,'Nucs')
        nuI = larmorfrq(Sys.Nucs,Exp.Field);
        for k = 1:numel(nuI)
          line([1 1]*abs(nuI(k)),ylim,'Color',[1 1 1]*0.8);
        end
        h = get(gca,'Children');
      end
      set(gca,'Children',h(end:-1:1));
    end

  else

    % Plot time-domain data matrix
    if plotFreqDomain
      subplot(1,3,1);
    end

    surf(x1,x2,real(info.td.'));
    view([0 90]);
    grid off
    shading flat
    axis tight
    box on
    set(gca,'Layer','top');
    title('Time domain (real part)');
    if ~isfield(Exp,'Dim1') || (strcmp(Exp.Dim1{1}(1),'d') && strcmp(Exp.Dim2{1}(1),'d'))
      axis equal
    end
    if ~isfield(Exp,'Dim1') || strcmp(Exp.Dim1{1}(1),'d')
      xlabel('{\itt}_1 (µs)');
    else
      xlabel(Exp.Dim1{1})
    end
    if ~isfield(Exp,'Dim2') || strcmp(Exp.Dim2{1}(1),'d')
      ylabel('{\itt}_2 (µs)');
    else
      xlabel(Exp.Dim2{1})
    end

    % Plot spectrum (first and second quadrant only)
    if plotFreqDomain
      subplot(1,3,[2 3]);

      if isfield(info,'f1') && isfield(info,'f2')
        fx1 = info.f1;
        fx2 = info.f2;
      else
        fx1 = fdaxis(Exp.dt(1),size(info.fd,1));
        if numel(Exp.dt)<2, Exp.dt(2) = Exp.dt(1); end
        fx2 = fdaxis(Exp.dt(2),size(info.fd,2));
      end
      fd = abs(info.fd);
      if isfield(Opt,'logplot') && Opt.logplot
        fd = log(fd);
        maxfd = max(max(fd));
        fd(fd<maxfd-6) = maxfd-6;
      end
      contour(fx1,fx2,fd.',20);
      shading flat
      axis equal tight
      grid on
      xti = xticks;
      set(gca,'YTick',xti);
      colorbar
      set(gca,'Layer','top');
      title('Frequency domain');
      xlabel('\nu_1 (MHz)');
      ylabel('\nu_2 (MHz)');

      fm1 = max(abs(fx1));
      fm2 = max(abs(fx2));
      fm = max(fm1,fm2);
      linecol = 'k';
      %line([-1 1]*fm,[0 0],'Color',linecol);
      line([0 0],[-1 1]*fm,'Color',linecol);
      line([-1 1]*fm,[-1 1]*fm,'Color',linecol,'LineStyle',':');
      line([1 -1]*fm,[-1 1]*fm,'Color',linecol,'LineStyle',':');
      ylim([0 1]*fm);
    end
  end

end

end
