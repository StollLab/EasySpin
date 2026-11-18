% saffron    Simulate pulse EPR signals and spectra
%
%     saffron(Sys,Exp)
%     saffron(Sys,Exp,Opt)
%     S = saffron(...)
%     [x,S] = saffron(...)
%     [x,S,info] = saffron(...)
%
%  Inputs:
%     Sys   ... spin system (electron spins and nuclei)
%     Exp   ... experimental parameters (time unit µs)
%     Opt   ... simulation options
%
%  Exp fields:
%     Field          magnetic field (mT)
%     mwFreq         microwave carrier frequency (GHz)
%     rfFreq         RF carrier frequency (MHz), for RF pulses
%     Sequence       pulse sequence: cell array of pulses (structures) and
%                    delays (µs); or the name of a predefined experiment,
%                    which sets the sequence and its incrementation scheme:
%                    '2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR',
%                    'DaviesENDOR' (macro parameters tau, T, t1, t2, dt,
%                    tp, Range, tprf)
%     nPoints, Dim1, Dim2, ...   incrementation scheme (indirect dimensions)
%     Pathways       coherence pathway selection: 'all' (default), 'echo'
%                    (echoes within the detection window), a string of
%                    electron coherence orders per pulse ('+0-', 'c' = + or
%                    -, '.' = any), or a numeric matrix (NaN = any)
%     PhaseCycle     phase cycle, one row per step: the phase of each pulse
%                    and the receiver phase (radians); the signal is the
%                    average over the steps
%     OnOff          pulses switched on and off (pulse numbers); the signal
%                    is the difference signal(on) - signal(off)
%     DetWindow      detection window (µs) relative to the end of the
%                    sequence: 0 (default, single point), [a b] (transient),
%                    or 'sequence' (whole sequence)
%     DetOperator    detection operators (sop strings or matrices), default
%                    S+ of all electrons
%     DetFreq, DetPhase, DetIntegrate, DetTimeStep
%     mwPolarization, ResonatorFrequency, ResonatorQL, FrequencyResponse,
%     ResonatorMode
%     SampleFrame, MolFrame, CrystalSymmetry, SampleRotation (powder by
%     default; Exp.MolFrame = [0 0 0] for a single orientation)
%
%  Pulse fields: Flip, tp (0 or absent: ideal pulse), Phase, Type,
%     Frequency (offset from the carrier, or [start end] for chirps),
%     Amplitude, Channel ('mw' or 'rf'), B1 (RF field, mT), and the
%     parameters of pulse(). For electron spins S > 1/2, microwave pulses
%     need tp and Amplitude (nutation frequency in MHz, pulse Hamiltonian
%     Amplitude*Sx); Flip, ideal pulses and predefined experiments are not
%     available.
%
%  Opt fields:
%     Method         'auto' (default), 'subspace', 'fullspace'
%     IdealPulses    treat finite pulses as ideal, with orientation
%                    selection by an excitation window of width 1/tp
%                    (subspace)
%     Relaxation     relaxation on/off, a single value or one per element
%                    of Exp.Sequence
%     GridSize, GridSymmetry, GridFrame, Sites
%     IntTimeStep, SimFreq (fullspace)
%     OffsetTolerance  error target of the average over the inhomogeneous
%                    line (Sys.lwAnisoG), default 1e-6 (subspace)
%     Binning        'auto' (default), true, false, or one true/false per
%                    dimension: delay dimensions computed by binning the
%                    peaks (fast, approximate) or point by point (exact)
%                    (subspace)
%     OriThreshold, Nuclei, Transitions, ProductRule, Expand (subspace)
%     separate, Window, ZeroFillFactor, Verbosity
%
%  Outputs:
%     x       axis of each indirect dimension (values of the first
%             parameter listed in Exp.DimN), and the transient time axis
%     S       simulated signal, dimensions [Dim1, Dim2, ..., transient
%             time, detection operator], singletons removed
%     info    structure with the expanded experiment (Exp), the engine
%             (Method, MethodReason), the processed spectrum (f, fd, or
%             f1, f2, fd), the start and end times of the elements of
%             Exp.Sequence on the time axis of x (EventTimes, fullspace),
%             the contributing coherence pathways (Pathways, in the
%             syntax of Exp.Pathways) and their echo times after the end
%             of the sequence (EchoTimes, µs, NaN without echo), the
%             number of resonance offsets in the average over the
%             inhomogeneous line (nOffsets, subspace), and the binned
%             dimensions (binnedDims, subspace)

function varargout = saffron(Sys,Exp,Opt)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker('internal'));

% Check Matlab version
warning(chkmlver);

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
logmsg(Opt.Verbosity);

if logmsg>=1, logmsg(1,'=begin=saffron====%s=================',char(datetime)); end
logmsg(2,'  log level %d',logmsg);
logmsg(1,'-general-----------------------------------------------');

% Process Opt.separate
if ~isfield(Opt,'separate'), Opt.separate = ''; end
[separateOutput,err] = parseoption(Opt,'separate',{'','components','transitions','orientations','sites'});
error(err);
if any(separateOutput==[3 4 5])
  error('saffron supports only Opt.separate = ''components''.');
end
separateComponentSpectra = separateOutput==2;

if ~isfield(Opt,'IsoCutoff'), Opt.IsoCutoff = 1e-4; end


%===============================================================================
% Experiment structure
%===============================================================================

DefaultExp.Temperature = [];
DefaultExp.Ordering = [];
DefaultExp.SampleRotation = [];
DefaultExp.SampleFrame = [];
DefaultExp.CrystalSymmetry = '';
DefaultExp.MolFrame = [];
Exp = adddefaults(Exp,DefaultExp);

% Removed and obsolete fields
if isfield(Exp,'Orientations')
  error('Exp.Orientations is no longer supported, use Exp.SampleFrame instead.');
end
if isfield(Exp,'CrystalOrientation')
  error('Exp.CrystalOrientation is no longer supported, use Exp.SampleFrame instead.');
end
obsoleteFields = {'Flip','Phase','tp','t','Inc'};
for index = 1:numel(obsoleteFields)
  if isfield(Exp,obsoleteFields{index}) && ~(ischar(Exp.Sequence) && strcmp(obsoleteFields{index},'tp'))
    error('Exp.%s is no longer supported, define pulses and delays in Exp.Sequence and increments in Exp.Dim1, Exp.Dim2 instead.',obsoleteFields{index});
  end
end
if isfield(Exp,'ExciteWidth')
  error(['Exp.ExciteWidth is no longer supported. Use finite pulses for orientation selection ' ...
    '(Exp.tp for predefined experiments, Pulse.tp otherwise), and Opt.IdealPulses = true for ' ...
    'the ideal-pulse approximation with an excitation window of width 1/tp.']);
end
if isfield(Exp,'Filter')
  error('Exp.Filter is no longer supported. Select coherence pathways with Exp.Pathways, e.g. ''+0-''.');
end
if isfield(Exp,'DetSequence')
  error(['Exp.DetSequence is no longer supported. Use Exp.DetWindow = ''sequence'' to detect ' ...
    'during the whole sequence; info.EventTimes gives the start and end of each element of ' ...
    'Exp.Sequence, for selecting the detected events.']);
end
if isfield(Exp,'lightBeam') && ~isempty(Exp.lightBeam)
  error('Photoselection (via Exp.lightBeam) is not supported.')
end
if isfield(Opt,'SimulationMode')
  error('Opt.SimulationMode is no longer supported. Use Opt.Method = ''auto'', ''subspace'' or ''fullspace''.');
end
if isfield(Opt,'lwOffset')
  error('Opt.lwOffset is no longer supported. The offset distribution is the inhomogeneous EPR line, Sys.lwAnisoG.');
end
if isfield(Opt,'EndorMethod')
  error('Opt.EndorMethod is no longer supported. ENDOR is simulated with the RF pulse in the sequence.');
end
if isfield(Opt,'nOffsets')
  error(['Opt.nOffsets is no longer supported. The number of resonance offsets is chosen ' ...
    'automatically; use Opt.OffsetTolerance to set the accuracy.']);
end
for f = {'TimeDomain','PartialIFFT','logplot','ExcOperator'}
  if isfield(Opt,f{1})
    error('Opt.%s is no longer supported.',f{1});
  end
end
if isfield(Opt,'DetOperator')
  error('Opt.DetOperator is not supported. Give the detection operators in Exp.DetOperator.');
end
if isfield(Opt,'StateTrajectories') && any(Opt.StateTrajectories)
  error(['Opt.StateTrajectories is not supported. To follow the spin dynamics, detect ' ...
    'expectation values over the whole sequence, e.g. Exp.DetOperator = {''x1'',''y1'',''z1''} ' ...
    'with Exp.DetWindow = ''sequence''.']);
end
if isfield(Exp,'Ordering') && ~isempty(Exp.Ordering)
  error('Partially ordered samples are not implemented in saffron.')
end

% Field
if ~isfield(Exp,'Field')
  error('Exp.Field is missing. Give a magnetic field in mT.');
end
if any(Exp.Field<0)
  error('Exp.Field cannot be negative.');
end

% Microwave carrier frequency, optional: empty is the same as not given
if isfield(Exp,'mwFreq') && isempty(Exp.mwFreq)
  Exp = rmfield(Exp,'mwFreq');
end

% Exp.Sequence is a required field
if ~isfield(Exp,'Sequence')
  error('Exp.Sequence is missing. Give a pulse sequence or the name of a predefined experiment.');
end

% Electron spins S > 1/2: with several EPR transitions, a flip angle has no
% unique meaning, so microwave pulses need their strength as Amplitude
if iscell(Sys), SysList = Sys; else, SysList = {Sys}; end
highSpin = false;
for c = 1:numel(SysList)
  highSpin = highSpin || (isfield(SysList{c},'S') && any(SysList{c}.S>1/2));
end
if highSpin
  checkHighSpinPulses(Exp,Opt);
end

% Predefined experiments are macros that expand into user-defined
% experiments (pulse sequence and incrementation scheme)
if ischar(Exp.Sequence) || isstring(Exp.Sequence)
  Exp.Sequence = char(Exp.Sequence);
  Exp = sf_predefinedexperiments(Exp);
  logmsg(1,'Predefined experiment %s, expanded.',Exp.Macro);
end

% Custom excitation operators are not supported (all pulses act on the
% total electron spin)
if iscell(Exp.Sequence)
  for k = find(cellfun(@isstruct,Exp.Sequence))
    if isfield(Exp.Sequence{k},'Operator')
      error('Pulse field Operator is not supported. Pulses act on all electron spins; select spins by their resonance frequencies instead.');
    end
  end
end

% RF pulses act on nuclei
if iscell(Exp.Sequence)
  hasRF = false;
  for k = find(cellfun(@isstruct,Exp.Sequence))
    for p = Exp.Sequence{k}(:).'
      hasRF = hasRF || (isfield(p,'Channel') && strcmpi(p.Channel,'rf'));
    end
  end
  if hasRF
    for c = 1:numel(SysList)
      if ~isfield(SysList{c},'Nucs') || isempty(SysList{c}.Nucs)
        error('The pulse sequence contains RF pulses, but the spin system has no nuclei (Sys.Nucs).');
      end
    end
  end
end

% Phase cycle and pulses switched on and off
if iscell(Exp.Sequence)
  Exp = parsePhaseCycle(Exp,sum(cellfun(@isstruct,Exp.Sequence)));
end

% Detection window, default: single point at the end of the sequence
if ~isfield(Exp,'DetWindow') || isempty(Exp.DetWindow)
  Exp.DetWindow = 0;
end
% Detection operators as a cell array; detection phases: one per operator
% (a single value applies to all operators). An empty DetOperator means
% the default (S+).
if isfield(Exp,'DetOperator') && ~isempty(Exp.DetOperator) && ~iscell(Exp.DetOperator)
  Exp.DetOperator = {Exp.DetOperator};
end
if isfield(Exp,'DetOperator') && ~isempty(Exp.DetOperator)
  nDetOps = numel(Exp.DetOperator);
else
  nDetOps = 1;
end
if isfield(Exp,'DetPhase') && ~isempty(Exp.DetPhase)
  if isscalar(Exp.DetPhase)
    Exp.DetPhase = repmat(Exp.DetPhase,1,nDetOps);
  elseif numel(Exp.DetPhase)~=nDetOps
    error('Exp.DetPhase must be a single value or have one value per detection operator in Exp.DetOperator (%d).',nDetOps);
  end
end

% Detection frame: one frequency for all detection operators
if isfield(Exp,'DetFreq') && ~isempty(Exp.DetFreq)
  if ~isnumeric(Exp.DetFreq) || ~isscalar(Exp.DetFreq) || ~isreal(Exp.DetFreq)
    error(['Exp.DetFreq must be a single frequency (in GHz). All signals and ' ...
      'density matrices are given in one detection frame; operators without ' ...
      'electron coherence, such as Sz, are not affected by it.']);
  end
  if Exp.DetFreq<0
    error('Exp.DetFreq cannot be negative.');
  end
end

% Temperature
if p_temperature(Exp)
  error('Exp.Temperature is not supported for pulse EPR simulations.');
end

% Detect sample type (disordered, partially ordered, crystal)
[Exp,Opt] = p_sampletype(Exp,Opt);
if Opt.partiallyOrderedSample
  error('Partially ordered samples are not implemented in saffron.')
end

%===============================================================================
% Options structure
%===============================================================================
if ~isfield(Opt,'GridSymmetry'), Opt.GridSymmetry = []; end
if ~isfield(Opt,'GridFrame'), Opt.GridFrame = []; end
if ~isfield(Opt,'Transitions'), Opt.Transitions = []; end
if ~isfield(Opt,'Sites'), Opt.Sites = []; end
if isfield(Opt,'nKnots')
  error('Options.nKnots is obsolete. Use Options.GridSize instead, e.g. Options.GridSize = 91.');
end
if isfield(Opt,'Relaxation') && ~isscalar(Opt.Relaxation) && numel(Opt.Relaxation)~=numel(Exp.Sequence)
  error('Opt.Relaxation must be a single value or have one entry per element of Exp.Sequence (%d).',numel(Exp.Sequence));
end
Opt.relaxationFromSys = ~isfield(Opt,'Relaxation');
if ~isfield(Opt,'GridSize'), Opt.GridSize = 30+1; end
if numel(Opt.GridSize)>1
  error('Only one number allowed in Opt.GridSize. saffron does not support interpolation.');
end
if Opt.GridSize<7
  error('Opt.GridSize must be at least 7. You gave %d.',Opt.GridSize);
end
if ~isfield(Opt,'Window'), Opt.Window = 'ham+'; end
if ~isfield(Opt,'ZeroFillFactor'), Opt.ZeroFillFactor = 2; end
if ~isfield(Opt,'IdealPulses'), Opt.IdealPulses = false; end
if ~isfield(Opt,'OffsetTolerance') || isempty(Opt.OffsetTolerance)
  Opt.OffsetTolerance = 1e-6;
end
if ~isnumeric(Opt.OffsetTolerance) || ~isscalar(Opt.OffsetTolerance) || ...
    ~isreal(Opt.OffsetTolerance) || ~(Opt.OffsetTolerance>0 && Opt.OffsetTolerance<1)
  error('Opt.OffsetTolerance must be a number between 0 and 1, e.g. 1e-6.');
end

% Choose the simulation engine
%===============================================================================
autoMethod = ~isfield(Opt,'Method') || isempty(Opt.Method) || strcmp(Opt.Method,'auto');
[Opt.Method,methodReason,seq,engineReasons] = sf_selectengine(Sys,Exp,Opt);
if autoMethod && ~isempty(engineReasons)
  logmsg(1,'Simulation engine: %s, since the subspace engine does not support',Opt.Method);
  for iReason = 1:numel(engineReasons)
    logmsg(1,'  - %s',engineReasons{iReason});
  end
else
  logmsg(1,'Simulation engine: %s %s',Opt.Method,methodReason);
end
subspaceEngine = strcmp(Opt.Method,'subspace');

if subspaceEngine
  Opt.fastSeq = seq;
  Opt.fastExp = Exp;
  pathwaySeq = seq;
else
  % Full-space engine: set up events and indirect dimensions
  ExpProp = Exp;
  if isfield(ExpProp,'DetOperator')
    Opt.DetOperator = ExpProp.DetOperator;
    ExpProp = rmfield(ExpProp,'DetOperator');
  end
  if (ischar(ExpProp.DetWindow) || isstring(ExpProp.DetWindow)) && ~strcmp(ExpProp.DetWindow,'sequence')
    error('Exp.DetWindow must be a number, two numbers, or ''sequence''.');
  end
  % Reference gyromagnetic ratio (MHz/mT) for nominal RF flip angles: the
  % first nucleus of the first component
  if iscell(Sys), Sys1 = Sys{1}; else, Sys1 = Sys; end
  if isfield(Sys1,'Nucs') && ~isempty(Sys1.Nucs)
    nuc1 = regexp(Sys1.Nucs,',','split');
    try
      Opt.rfGammaRef = abs(nucgval(strtrim(nuc1{1})))*nmagn/1e3/planck/1e6;
    catch
    end
  end
  % Propagation runs of the phase cycle, and events and indirect dimensions
  % with the wave forms of all settings of each pulse in the phase cycle
  [cycleSettings,cycleRuns] = phaseCycleRuns(ExpProp,sum(cellfun(@isstruct,Exp.Sequence)));
  logmsg(1,'  phase cycle: %d propagation run(s)',numel(cycleRuns));
  ExpProp = rmfield(ExpProp,intersect(fieldnames(ExpProp),{'PhaseCycle','OnOff'}));
  ExpProp.pulseCycles = cellfun(@(s)[s(:,1) ones(size(s,1),1) s(:,2)],cycleSettings,'UniformOutput',false);
  [Events,Vary,Opt] = sf_buildevents(ExpProp,Opt);
  Opt.propEvents = Events;
  Opt.propVary = Vary;
  Opt.propCycleRuns = cycleRuns;
  % Pathways of a single EPR transition (info.Pathways): one electron spin
  % with S = 1/2 in each component
  if iscell(Sys), SysList = Sys; else, SysList = {Sys}; end
  spinHalf = all(cellfun(@(S_) ~isfield(S_,'S') || isequal(S_.S,1/2),SysList));
  [Opt.propPathways,pathwaySeq] = sf_pathwayorders(Exp,Opt.IdealPulses,spinHalf);
end

%===============================================================================
% Loop over components and isotopologues
%===============================================================================
autoRange = false;
includeInverseDomain = true;
[~,data,info] = compisoloop(@saffron_single,Sys,Exp,Opt,autoRange,true,separateComponentSpectra,includeInverseDomain);

%===============================================================================
% Axes, output arrangement, processing
%===============================================================================
x = dimensionAxes(Exp);
if ~subspaceEngine && isfield(info,'transientAxis') && ~isempty(info.transientAxis)
  x{end+1} = info.transientAxis;
end
if isscalar(x), x = x{1}; end
if isempty(x), x = []; end

info.Exp = Exp;
info.Method = Opt.Method;
info.MethodReason = methodReason;

% Contributing coherence pathways (alpha and beta populations both have
% order 0), in the order of their echo times (pathways without echo last)
info.Pathways = [];
info.EchoTimes = [];
if ~isempty(pathwaySeq)
  contributing = pathwaySeq.pathwayprefactor~=0;
  [co,iFirst] = unique(pathwaySeq.co(contributing,:),'rows','stable');
  echoTimes = pathwaySeq.echoTimes(contributing);
  [info.EchoTimes,idx] = sort(echoTimes(iFirst));
  symbols = '-0+';
  info.Pathways = reshape(symbols(co(idx,:)+2),size(co));
end
if isfield(info,'transientAxis'), info = rmfield(info,'transientAxis'); end
if isfield(info,'decayRates')
  decayAdded = ~isempty(info.decayRates) && any(info.decayRates(:)>0);
else
  decayAdded = false;
end
if ~subspaceEngine
  decayAdded = isfield(info,'relaxation') && info.relaxation;
end
info = processSignal(x,data,Exp,Opt,decayAdded,info);

switch nargout
  case 0
    saffron_plot(x,data,info,Exp);
  case 1
    varargout = {data};
  case 2
    varargout = {x,data};
  case 3
    varargout = {x,data,info};
end

%===============================================================
% Report performance
%===============================================================
endTime = datetime;
elapsedtime = endTime-startTime;
elapsedtime.Format = 'hh:mm:ss.SSS';
logmsg(1,'saffron took %s',elapsedtime);

if logmsg>=1, logmsg(1,'=end=saffron======%s=================\n',datetime); end

end


%===============================================================================
% saffron_single  Simulate a single isotopologue of a single component
%
%   [x_out,y_out,info] = saffron_single(Sys,Exp,Opt)
%
%   Validates the spin system, sets up the orientations, and runs the
%   simulation engine. Called by saffron via compisoloop.
%
%   Input:
%     Sys    spin system of a single isotopologue
%     Exp    experiment, parsed by saffron
%     Opt    options, parsed by saffron, including
%       Method             'subspace' or 'fullspace'
%       relaxationFromSys  true if the user did not give Opt.Relaxation
%       fastSeq, fastExp   parsed sequence and experiment (subspace engine)
%       propEvents, propVary, propCycleRuns, propPathways
%                          events and indirect dimensions with all
%                          settings of the pulses in the phase cycle, the
%                          propagation runs of the phase cycle (see
%                          phaseCycleRuns), and electron coherence orders
%                          of the selected pathways (full-space engine)
%
%   Output:
%     x_out  always empty (saffron computes the axes)
%     y_out  signal, see sf_engine_subspace and sf_engine_fullspace
%     info   structure from the engine, with additional fields
%       SinglePointDetection  always true
%       relaxation            true if relaxation is on (full-space engine)
function [x_out,y_out,info] = saffron_single(Sys,Exp,Opt)

subspaceEngine = strcmp(Opt.Method,'subspace');
x_out = [];

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
if ~isempty(Sys.StrainData.Q)
  error('saffron does not support strains (Sys.StrainPars).');
end
if isfield(Sys,'T1T2')
  error('Sys.T1T2 is not supported. Provide relaxation times in Sys.T1 and Sys.T2.');
end
if isfield(Sys,'n') && any(Sys.n~=1)
  error('saffron does not support equivalent nuclei (Sys.n). Specify equivalent nuclei explicitly.');
end

% Relaxation: on if not given and the spin system has nonzero relaxation times
relaxation = (isfield(Sys,'T1') && any(Sys.T1(:))) || (isfield(Sys,'T2') && any(Sys.T2(:)));
if Opt.relaxationFromSys
  Opt.Relaxation = relaxation;
end

% Set up grid
if ~isfield(Exp,'OriWeights')
  [Exp,Opt] = p_gridsetup(Sys,Exp,Opt);
end

% Process crystal orientations, crystal symmetry, and frame transforms
[Orientations,nOrientations,nSites,averageOverChi] = p_crystalorientations(Exp,Opt);
if numel(Exp.OriWeights)~=nOrientations
  Exp.OriWeights = repmat(Exp.OriWeights,1,nSites);
end
% Total orientation weight 1
Exp.OriWeights = Exp.OriWeights/sum(Exp.OriWeights);

if subspaceEngine
  [y_out,info] = sf_engine_subspace(Sys,Exp,Opt,Orientations);
else
  [y_out,info] = sf_engine_fullspace(Sys,Exp,Opt,Orientations,nOrientations,averageOverChi);
  info.relaxation = any(Opt.Relaxation);
end

% ENDOR line broadening along RF carrier sweeps
rfDims = rfFrequencyDims(Exp);
if isscalar(rfDims) && isvector(y_out) && isfield(Sys,'lwEndor') && any(Sys.lwEndor>0)
  df = Exp.(sprintf('Dim%d',rfDims)){1,2};
  if ~isscalar(df), df = df(2)-df(1); end
  y_out = convspec(y_out,abs(df),Sys.lwEndor(1));
end

info.SinglePointDetection = true;

end


%===============================================================================
% sf_pathwayorders  Electron coherence orders of the selected pathways
%
%   orders = sf_pathwayorders(Exp,idealPulses)
%
%   Converts the pathway selection in Exp.Pathways to the numeric form used
%   by the full-space engine.
%
%   Input:
%     Exp     experiment; uses
%       Pathways  'all', 'echo', a string of electron coherence orders per
%                 pulse ('+', '-', '0', 'c' = + or -, '.' = any), or a
%                 numeric matrix
%       Sequence  pulse sequence
%     idealPulses passed on to sf_parsesequence
%     spinHalf    true if the spin system has one electron spin with S = 1/2
%
%   Output:
%     orders  electron coherence order after each pulse (+1, -1, 0, NaN =
%             any, Inf = +1 or -1), one row per pathway and one column per
%             pulse; empty if there is no selection
%     seq     parsed sequence with the pathways of the single EPR
%             transition (from sf_parsesequence), for info.Pathways; empty
%             if the spin system has other electron spins or the selection
%             has coherence orders beyond 1
function [orders,seq] = sf_pathwayorders(Exp,idealPulses,spinHalf)
orders = [];
seq = [];
pw = 'all';
if isfield(Exp,'Pathways') && ~isempty(Exp.Pathways), pw = Exp.Pathways; end
if isstring(pw), pw = char(pw); end

% sf_parsesequence checks the selection in the same way for both engines,
% and gives the pathways of a single EPR transition. It only knows coherence
% orders up to 1, higher ones need S>1/2 or several spins.
higherOrders = isnumeric(pw) && any(abs(pw(:))>1 & isfinite(pw(:)));
if higherOrders
  nPulses = sum(cellfun(@isstruct,Exp.Sequence));
  if size(pw,2)~=nPulses
    error('Exp.Pathways has %d entries per pathway, but the sequence has %d pulses.',size(pw,2),nPulses);
  end
elseif ~strcmp(pw,'all') || (spinHalf && any(cellfun(@isstruct,Exp.Sequence)))
  seq = sf_parsesequence(Exp,idealPulses);
end

if strcmp(pw,'all')
  % no selection
elseif strcmp(pw,'echo')
  orders = unique(seq.echoOrders,'rows');
elseif ischar(pw)
  orders = nan(size(pw));
  orders(pw=='+') = +1;
  orders(pw=='-') = -1;
  orders(pw=='0') = 0;
  orders(pw=='c') = Inf;
else
  orders = pw;
end
if ~spinHalf, seq = []; end
end

%===============================================================================
% parsePhaseCycle  Check the phase cycle and the pulses switched on and off
%
%   Exp = parsePhaseCycle(Exp,nPulses)
%
%   Checks Exp.PhaseCycle and Exp.OnOff.
%
%   Input:
%     Exp      experiment; uses PhaseCycle, OnOff
%     nPulses  number of pulses in Exp.Sequence
%
%   Output:
%     Exp      experiment with OnOff as a row of pulse numbers, if given
function Exp = parsePhaseCycle(Exp,nPulses)

if isfield(Exp,'PhaseCycle') && ~isempty(Exp.PhaseCycle)
  PC = Exp.PhaseCycle;
  if iscell(PC)
    error(['Exp.PhaseCycle with one phase cycle per pulse (cell array) is no longer supported. ' ...
      'Give the phase cycle as a matrix with one row per step, containing the phases of all ' ...
      'pulses and the receiver phase (in radians).']);
  end
  if ~isnumeric(PC) || ~isreal(PC) || ~ismatrix(PC) || any(~isfinite(PC(:)))
    error('Exp.PhaseCycle must be a matrix of phases (in radians).');
  end
  if size(PC,2)~=nPulses+1
    error(['Exp.PhaseCycle must have one row per step and %d columns: the phases of the ' ...
      '%d pulses and the receiver phase (in radians). It has %d columns.'],...
      nPulses+1,nPulses,size(PC,2));
  end
end

if isfield(Exp,'OnOff') && ~isempty(Exp.OnOff)
  onOff = Exp.OnOff;
  if ~isnumeric(onOff) || ~isvector(onOff) || any(onOff<1 | onOff>nPulses | mod(onOff,1)~=0)
    error('Exp.OnOff must contain pulse numbers between 1 and %d.',nPulses);
  end
  if numel(unique(onOff))<numel(onOff)
    error('Exp.OnOff contains a pulse number more than once.');
  end
  Exp.OnOff = onOff(:).';
end

end

%===============================================================================
% phaseCycleRuns  Propagation runs of the phase cycle (full-space engine)
%
%   [settings,runs] = phaseCycleRuns(Exp,nPulses)
%
%   The signal is a weighted sum over acquisition steps: the steps of the
%   phase cycle, each with the weight exp(-1i*receiver phase)/nSteps, and,
%   with Exp.OnOff, the same steps with the listed pulses switched off and
%   the opposite weight. Steps with identical pulse settings are combined.
%   Phases are compared modulo 2*pi, and the phase of a pulse that is
%   switched off is ignored.
%
%   If the steps form a nested phase cycle, i.e. all combinations of the
%   settings (phase, on/off) of the individual pulses, with weights that are
%   products of one factor per setting of each pulse, the sum is computed in
%   a single run: since propagation is linear, sf_propagate adds up the
%   weighted states of the settings of each pulse right after the pulse.
%   Otherwise, each step with nonzero weight is propagated in a run of its
%   own.
%
%   Input:
%     Exp      experiment; uses PhaseCycle and OnOff (from parsePhaseCycle)
%     nPulses  number of pulses in Exp.Sequence
%
%   Output:
%     settings  cell array with one element per pulse: its settings in the
%               phase cycle, rows [phase scale] (phase in radians, added to
%               Pulse.Phase; scale 1 for on and 0 for off)
%     runs      structure array, one element per run, with fields
%       rows      cell array with one element per pulse: the settings used
%                 in the run (row indices into settings{p})
%       factors   cell array with one element per pulse: the weights of
%                 these settings (column)
%       weight    weight of the signal of the run
function [settings,runs] = phaseCycleRuns(Exp,nPulses)

% Acquisition steps: phases and scales of the pulses, signal weights
if isfield(Exp,'PhaseCycle') && ~isempty(Exp.PhaseCycle)
  phases = Exp.PhaseCycle(:,1:nPulses);
  weights = exp(-1i*Exp.PhaseCycle(:,end))/size(Exp.PhaseCycle,1);
else
  phases = zeros(1,nPulses);
  weights = 1;
end
scales = ones(size(phases));
if isfield(Exp,'OnOff') && ~isempty(Exp.OnOff)
  off = scales;
  off(:,Exp.OnOff) = 0;
  phases = [phases; phases];
  scales = [scales; off];
  weights = [weights; -weights];
end

% No pulses: a single run
if nPulses==0
  settings = cell(1,0);
  runs = struct('rows',{cell(1,0)},'factors',{cell(1,0)},'weight',sum(weights));
  return
end

% Phases modulo 2*pi; the phase of a pulse that is switched off has no
% effect
tol = 1e-10;
phases = mod(phases,2*pi);
phases(abs(phases-2*pi)<tol | scales==0) = 0;

% Settings of each pulse, and the setting of each pulse in each step;
% identical steps are combined, and weights that cancel up to rounding
% errors are set to zero
settings = cell(1,nPulses);
idx = zeros(numel(weights),nPulses);
for p = 1:nPulses
  [settings{p},~,idx(:,p)] = uniquetol([phases(:,p) scales(:,p)],tol,'ByRows',true,'DataScale',1);
end
weightTol = 1e-12*sum(abs(weights));
[idx,~,iStep] = unique(idx,'rows');
weights = accumarray(iStep,weights);
weights(abs(weights)<=weightTol) = 0;
nSettings = cellfun(@(s)size(s,1),settings);

% Nested phase cycle: the array of the weights of all combinations of the
% settings (one dimension per pulse) is an outer product of one factor
% vector per pulse
nested = size(idx,1)==prod(nSettings);
if nested
  W = accumarray([idx ones(size(idx,1),1)],weights,[nSettings 1]);
  [~,iRef] = max(abs(W(:)));
  ref = cell(1,nPulses);
  [ref{:}] = ind2sub([nSettings 1],iRef);
  factors = cell(1,nPulses);
  Wnested = W(iRef);
  for p = 1:nPulses
    sub = ref;
    sub{p} = ':';
    factors{p} = reshape(W(sub{:}),[],1)/W(iRef);
    Wnested = Wnested.*reshape(factors{p},[ones(1,p-1) nSettings(p) 1]);
  end
  nested = W(iRef)~=0 && norm(Wnested(:)-W(:))<=1e-10*norm(W(:));
end

if nested
  % settings with factor zero do not contribute (the setting at iRef has
  % factor 1, so each pulse keeps at least one)
  rows = cell(1,nPulses);
  for p = 1:nPulses
    rows{p} = find(factors{p}~=0);
    factors{p} = factors{p}(rows{p});
  end
  runs = struct('rows',{rows},'factors',{factors},'weight',W(iRef));
else
  keep = find(weights~=0);
  if isempty(keep), keep = 1; end  % no signal: one run with zero weight
  runs = struct('rows',cell(1,numel(keep)),'factors',[],'weight',[]);
  for r = 1:numel(keep)
    runs(r).rows = num2cell(idx(keep(r),:));
    runs(r).factors = num2cell(ones(1,nPulses));
    runs(r).weight = weights(keep(r));
  end
end

% Keep only the settings that are used in the runs
for p = 1:nPulses
  used = unique(cell2mat(arrayfun(@(r)r.rows{p}(:),runs(:),'UniformOutput',false)));
  newRow = zeros(nSettings(p),1);
  newRow(used) = 1:numel(used);
  settings{p} = settings{p}(used,:);
  for r = 1:numel(runs)
    runs(r).rows{p} = newRow(runs(r).rows{p});
  end
end

end

%===============================================================================
% rfFrequencyDims  Find the indirect dimensions that sweep the RF carrier
%
%   d = rfFrequencyDims(Exp)
%
%   Input:
%     Exp  experiment; uses nPoints, Dim1, Dim2, ...
%
%   Output:
%     d    indices of the dimensions whose first target is the RF carrier
%          frequency ('rfFreq'); empty if there are none
function d = rfFrequencyDims(Exp)
d = [];
if ~isfield(Exp,'nPoints'), return; end
for k = 1:numel(Exp.nPoints)
  f = sprintf('Dim%d',k);
  if isfield(Exp,f) && strcmp(strtrim(Exp.(f){1,1}),'rfFreq')
    d(end+1) = k; %#ok<AGROW>
  end
end
end

%===============================================================================
% dimensionAxes  Compute the axes of the indirect dimensions
%
%   x = dimensionAxes(Exp)
%
%   Input:
%     Exp  experiment; uses nPoints, Dim1, Dim2, ..., Sequence, and the
%          fields rfFreq, mwFreq, Field if they are incremented
%
%   Output:
%     x    cell array with one axis per indirect dimension: the values of
%          the first parameter listed in each Exp.DimN (its initial value
%          plus the increments); empty if there are no indirect dimensions
function x = dimensionAxes(Exp)
x = {};
if ~isfield(Exp,'nPoints') || isempty(Exp.nPoints), return; end
pulseIdx = find(cellfun(@isstruct,Exp.Sequence));
delayIdx = find(~cellfun(@isstruct,Exp.Sequence));
for d = 1:numel(Exp.nPoints)
  Dim = Exp.(sprintf('Dim%d',d));
  target = regexp(Dim{1,1},',','split');
  target = strtrim(target{1});
  v = Dim{1,2};
  n = Exp.nPoints(d);
  if isscalar(v), vals = v*(0:n-1); else, vals = v(:).'; end
  x0 = 0;
  if ~isempty(regexp(target,'^d\d+$','once'))
    x0 = Exp.Sequence{delayIdx(str2double(target(2:end)))};
  elseif any(strcmp(target,{'rfFreq','mwFreq','Field'}))
    x0 = Exp.(target);
  else
    tok = regexp(target,'^p(\d+)\.(\w+)','tokens','once');
    if ~isempty(tok)
      p_ = Exp.Sequence{pulseIdx(str2double(tok{1}))};
      if isfield(p_,tok{2}) && isnumeric(p_.(tok{2})) && ~isempty(p_.(tok{2}))
        x0 = p_.(tok{2})(1);
      end
    end
  end
  x{d} = x0 + vals; %#ok<AGROW>
end
end

%===============================================================================
% processSignal  Compute the spectrum of a simulated signal
%
%   info = processSignal(x,data,Exp,Opt,decayAdded,info)
%
%   Applies baseline correction, apodization and Fourier transform to
%   single-point signals along one or two delay dimensions. For 1D RF
%   frequency sweeps (ENDOR), the signal is returned as the spectrum. For
%   all other signals, info is returned unchanged.
%
%   Input:
%     x           axis, or cell array of axes, see dimensionAxes
%     data        simulated signal
%     Exp         experiment; uses nPoints, Dim1, Dim2
%     Opt         options; uses Window, ZeroFillFactor
%     decayAdded  true if the signal contains relaxation decays, which
%                 are then removed by the baseline correction
%     info        structure to add the results to
%
%   Output:
%     info        input structure, with additional fields
%       f         frequency axis (MHz) of a 1D spectrum, RF axis for ENDOR
%       f1, f2    frequency axes (MHz) of a 2D spectrum
%       td        time-domain signal (not for ENDOR)
%       fd        spectrum
function info = processSignal(x,data,Exp,Opt,decayAdded,info)
if isempty(x) || isscalar(data) || ~isfield(Exp,'nPoints') || isempty(Exp.nPoints), return; end
if ~iscell(x), x = {x}; end
nDims = numel(Exp.nPoints);
if ~isempty(rfFrequencyDims(Exp)) && nDims==1
  info.f = x{1};
  info.fd = data;
  return
end
if numel(x)~=nDims || nDims>2, return; end
% all dimensions must be delay increments
for d = 1:nDims
  Dim = Exp.(sprintf('Dim%d',d));
  if any(cellfun(@isempty,regexp(Dim(:,1),'^(d\d+)(,d\d+)*$','once'))) || ~isscalar(Dim{1,2})
    return
  end
end
if size(data,1)>1 && nDims==1, return; end  % separate components
dt = zeros(1,nDims);
for d = 1:nDims
  dt(d) = Exp.(sprintf('Dim%d',d)){1,2};
end
td = data;
if nDims==1
  if decayAdded && numel(td)>3
    [~,~,tdfit] = exponfit(x{1},td,1);
    tdx = td - tdfit;
  else
    tdx = td;
  end
  tdx = tdx - mean(tdx);
  win = apowin(Opt.Window,numel(tdx)).';
  tdx = tdx.*win;
  fd = fftshift(fft(tdx,Opt.ZeroFillFactor*numel(tdx)));
  info.f = fdaxis(dt(1),numel(fd));
else
  if decayAdded, order = 2; else, order = 0; end
  tdx = basecorr(td,1,order);
  tdx = basecorr(tdx,2,order);
  w1 = apowin(Opt.Window,size(td,1));
  w2 = apowin(Opt.Window,size(td,2));
  fd = fftshift(fftn(tdx.*(w1*w2.'),Opt.ZeroFillFactor*size(td)));
  info.f1 = fdaxis(dt(1),size(fd,1));
  info.f2 = fdaxis(dt(2),size(fd,2));
end
if max(abs(fd(:)))<1e-300, fd = fd*0; end
info.td = td;
info.fd = fd;
end

%===============================================================================
% saffron_plot  Plot the simulated signal and spectrum
%
%   saffron_plot(x,data,info,Exp)
%
%   Input:
%     x     axis, or cell array of axes: indirect dimensions, then the
%           transient time axis (if any)
%     data  simulated signal, as returned by saffron
%     info  structure with the processed spectrum (f or f1, f2, and fd),
%           see processSignal
%     Exp   experiment; uses Macro, DetOperator, Field, nPoints, Dim1, ...
%
%   Output:
%     none; plots into the current figure. Single points are printed, and
%     outputs with more than two axes are not plotted.
function saffron_plot(x,data,info,Exp)

logmsg(1,'Plotting...');
if isfield(Exp,'Macro'), name = Exp.Macro; else, name = 'User-defined experiment'; end

% Axes: indirect dimensions, then the transient time axis (if any)
if isempty(x), ax = {}; elseif iscell(x), ax = x; else, ax = {x}; end
nAx = numel(ax);
xl = axisLabels(Exp,nAx);
if isfield(Exp,'DetOperator') && ~isempty(Exp.DetOperator)
  ops = Exp.DetOperator;
  yLabel = 'expectation value';
else
  ops = {};
  yLabel = 'echo amplitude';
end
nDet = max(numel(ops),1);

% Single points and outputs with more than two axes are not plotted
if nAx==0
  fprintf('%s: signal = %s\n',name,mat2str(data,4));
  return
elseif nAx>2
  fprintf('%s: no plot for outputs with more than two axes.\n',name);
  return
end
clf

% Several detection operators: one subplot per operator
if nDet>1
  sz = cellfun(@numel,ax);
  D = reshape(data,[sz nDet]);
  for k = 1:nDet
    subplot(nDet,1,k);
    if ischar(ops{k}), opName = ['\langle' ops{k} '\rangle']; else, opName = sprintf('operator %d',k); end
    if nAx==1
      d = D(:,k).';
      if isreal(d) || max(abs(imag(d)))<1e-6*max(abs(d))
        plot(ax{1},real(d));
      else
        plot(ax{1},real(d),ax{1},imag(d));
        legend('real','imag'); legend boxoff
      end
      axis tight
      xlabel(xl{1}); ylabel(yLabel);
      title(opName);
    else
      surf(ax{1},ax{2},real(D(:,:,k)).');
      view([0 90]); shading flat; axis tight; box on
      set(gca,'Layer','top');
      xlabel(xl{1}); ylabel(xl{2});
      title([opName ' (real part)']);
    end
  end
  return
end

isENDOR = ~isempty(rfFrequencyDims(Exp));
if isENDOR && isvector(data)
  plot(x,real(data));
  xlabel('RF frequency (MHz)');
  ylabel('ENDOR effect');
  title(sprintf('%s, %g mT',name,Exp.Field));
  axis tight
  return
end

if iscell(x) && numel(x)==2 && ismatrix(data) && min(size(data))>1
  % 2D
  hasFD = isfield(info,'fd') && isfield(info,'f1');
  if hasFD, subplot(1,3,1); end
  surf(x{1},x{2},real(data.'));
  view([0 90]); shading flat; axis tight; box on
  set(gca,'Layer','top');
  xlabel(xl{1}); ylabel(xl{2});
  title('Time domain (real part)');
  if hasFD
    subplot(1,3,[2 3]);
    contour(info.f1,info.f2,abs(info.fd).',20);
    axis equal tight; grid on; colorbar
    set(gca,'Layer','top');
    fm = max(max(abs(info.f1)),max(abs(info.f2)));
    line([0 0],[-1 1]*fm,'Color','k');
    line([-1 1]*fm,[-1 1]*fm,'Color','k','LineStyle',':');
    line([1 -1]*fm,[-1 1]*fm,'Color','k','LineStyle',':');
    ylim([0 1]*fm);
    xlabel('\nu_1 (MHz)'); ylabel('\nu_2 (MHz)');
    title(sprintf('%s, spectrum',name));
  end
  return
end

if iscell(x), x = x{1}; end
if ~isvector(data)
  plot(x,real(data.'));
  xlabel(xl{1});
  ylabel(yLabel);
  title(name);
  return
end
hasFD = isfield(info,'fd') && isfield(info,'f');
if hasFD, subplot(2,1,1); end
if isreal(data) || max(abs(imag(data)))<1e-6*max(abs(data))
  plot(x,real(data));
else
  plot(x,real(data),x,imag(data));
  legend('real','imag'); legend boxoff
end
axis tight
xlabel(xl{1});
ylabel(yLabel);
title(sprintf('%s, time domain',name));
if hasFD
  subplot(2,1,2);
  idx = find(info.f>=0);
  plot(info.f(idx),abs(info.fd(idx)),info.f(idx),real(info.fd(idx)));
  legend('abs','real'); legend boxoff
  axis tight
  xlabel('\nu (MHz)');
  ylabel('intensity (arb.u.)');
  title('Spectrum');
end
end

% axisLabels  Axis labels for plotting
%
%   xl = axisLabels(Exp,nAx)
%
%   Input:
%     Exp  experiment; uses nPoints, Dim1, Dim2, ...
%     nAx  number of axes
%
%   Output:
%     xl   cell array of label strings, at least two: indirect dimensions,
%          then the transient time
function xl = axisLabels(Exp,nAx)
xl = repmat({'t (µs)'},1,max(nAx,2));
if ~isfield(Exp,'nPoints'), return; end
for d = 1:min(numel(Exp.nPoints),nAx)
  Dim = Exp.(sprintf('Dim%d',d));
  target = strtrim(Dim{1,1});
  if strcmp(target,'rfFreq')
    xl{d} = 'RF frequency (MHz)';
  elseif target(1)=='d'
    xl{d} = sprintf('%s (µs)',target);
  else
    xl{d} = target;
  end
end
end

%===============================================================================
% checkHighSpinPulses  Check microwave pulses for electron spins S > 1/2
%
%   checkHighSpinPulses(Exp,Opt)
%
%   The pulse Hamiltonian is H1 = nu1*Sx, with nu1 = Amplitude (MHz), and
%   transition a<->b nutates with nu1*|<a|S+|b>|. A flip angle (Flip, or
%   Qcrit, which sets the amplitude by a spin-1/2 criterion) has no unique
%   meaning with several transitions, so only Amplitude is accepted, and
%   only for finite pulses.
%
%   Input:
%     Exp  experiment; uses Sequence and Dim1, Dim2, ...
%     Opt  options; uses IdealPulses
%
%   Output:
%     none; throws an error if a microwave pulse does not comply
function checkHighSpinPulses(Exp,Opt)

if ischar(Exp.Sequence) || isstring(Exp.Sequence)
  error(['Predefined experiments are not available for electron spins S > 1/2. ' ...
    'Give the pulse sequence explicitly in Exp.Sequence, with the pulse strength ' ...
    'as Amplitude (nutation frequency in MHz).']);
end
if isfield(Opt,'IdealPulses') && Opt.IdealPulses
  error('Opt.IdealPulses is not supported for electron spins S > 1/2. Use finite pulses.');
end
msgAmplitude = ['give its strength as Amplitude (nutation frequency nu1 in MHz, ' ...
  'pulse Hamiltonian nu1*Sx).'];

if ~iscell(Exp.Sequence), return; end  % checked in the engines
pulseIdx = find(cellfun(@isstruct,Exp.Sequence));
isRF = false(1,numel(pulseIdx));
for k = 1:numel(pulseIdx)
  pulses = Exp.Sequence{pulseIdx(k)};
  for p = pulses(:).'
    if isfield(p,'Channel') && strcmpi(p.Channel,'rf')
      isRF(k) = true;
      continue
    end
    if isfield(p,'IQ'), continue; end  % user-defined waveform, in MHz
    for f = {'Flip','Qcrit'}
      if isfield(p,f{1}) && ~isempty(p.(f{1}))
        error(['Pulse %d: %s is not accepted for electron spins S > 1/2, since a flip angle ' ...
          'has no unique meaning with several EPR transitions. Instead, %s'],k,f{1},msgAmplitude);
      end
    end
    if ~isfield(p,'Amplitude') || isempty(p.Amplitude)
      error('Pulse %d: for electron spins S > 1/2, %s',k,msgAmplitude);
    end
    if ~isfield(p,'tp') || isempty(p.tp) || all(p.tp==0)
      error('Pulse %d: ideal pulses (no tp) are not supported for electron spins S > 1/2. Give a pulse length tp and an Amplitude.',k);
    end
  end
end

% Indirect dimensions over flip angles of microwave pulses
dimFields = fieldnames(Exp);
dimFields = dimFields(~cellfun(@isempty,regexp(dimFields,'^Dim\d+$','once')));
for d = 1:numel(dimFields)
  Dim = Exp.(dimFields{d});
  if ~iscell(Dim), continue; end
  for iLine = 1:size(Dim,1)
    targets = strtrim(regexp(Dim{iLine,1},',','split'));
    for t = 1:numel(targets)
      tok = regexp(targets{t},'^p(\d+)\.(Flip|Qcrit)$','tokens','once');
      if isempty(tok), continue; end
      iPulse = str2double(tok{1});
      if iPulse>=1 && iPulse<=numel(isRF) && isRF(iPulse), continue; end
      error('Exp.%s: for electron spins S > 1/2, increment p%d.Amplitude instead of p%d.%s.',...
        dimFields{d},iPulse,iPulse,tok{2});
    end
  end
end

end
