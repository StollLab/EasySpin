% saffron    Simulate pulse EPR signals and spectra
%
%     saffron(Sys,Exp)
%     saffron(Sys,Exp,Opt)
%     S = saffron(...)
%     [x,S] = saffron(...)
%     [x,S,info] = saffron(...)
%
%  Inputs:
%     Sys   ... spin system with electron spin and nuclei
%     Exp   ... experimental parameters (time unit µs)
%     Opt   ... simulation options
%
%  Exp fields:
%     Field          magnetic field (mT)
%     mwFreq         microwave carrier frequency (GHz)
%     rfFreq         RF carrier frequency (MHz), for RF pulses
%     Sequence       pulse sequence: cell array of pulses (structures) and
%                    delays (µs), or the name of a predefined experiment:
%                    '2pESEEM','3pESEEM','4pESEEM','HYSCORE','MimsENDOR',
%                    'DaviesENDOR' (macro parameters tau, T, t1, t2, dt,
%                    tp, Range, tprf)
%     nPoints, Dim1, Dim2, ...   indirect dimensions
%     Pathways       coherence pathway selection: 'all' (default), 'echo',
%                    a string of electron coherence orders per pulse
%                    ('+0-', '.' = any), or a numeric matrix (NaN = any)
%     PhaseCycle     phase cycle per pulse, rows [phase weight] or
%                    [phase weight scale]
%     DetWindow      detection window (µs) relative to the end of the
%                    sequence: 0 (default, single point), [a b] (transient),
%                    or 'sequence' (whole sequence)
%     DetOperator, DetFreq, DetPhase, DetIntegrate, DetTimeStep
%     mwPolarization, ResonatorFrequency, ResonatorQL, FrequencyResponse
%     SampleFrame, MolFrame, CrystalSymmetry, SampleRotation
%
%  Pulse fields: Flip, tp (0 or absent: ideal pulse), Phase, Type,
%     Frequency (offset from the carrier, or [start end] for chirps),
%     Amplitude, Channel ('mw' or 'rf'), B1 (RF field, mT), and the
%     parameters of pulse().
%
%  Opt fields:
%     Method         'auto' (default), 'subspace', 'fullspace'
%     IdealPulses    treat finite pulses as ideal, with orientation
%                    selection by an excitation window of width 1/tp
%     GridSize, GridSymmetry, GridFrame, nOffsets
%     IntTimeStep, SimFreq (fullspace)
%     OriThreshold, Nuclei, Transitions, ProductRule (subspace)
%     Relaxation, separate, Window, ZeroFillFactor, Verbosity
%
%  Outputs:
%     x       axis of each indirect dimension (values of the first
%             parameter listed in Exp.DimN), and the transient time axis
%     S       simulated signal, dimensions [Dim1, Dim2, ..., transient
%             time, detection operator], singletons removed
%     info    structure with the expanded experiment (Exp), the engine
%             (Method, MethodReason), the processed spectrum (f, fd), etc.

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
for f = {'TimeDomain','PartialIFFT','logplot','ExcOperator'}
  if isfield(Opt,f{1})
    error('Opt.%s is no longer supported.',f{1});
  end
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
  error('Exp.Sequence is missing. Provide a predefined or user-defined pulse sequence.');
end

% Predefined experiments are macros that expand into user-defined sequences
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
if isfield(Opt,'Relaxation') && ~isscalar(Opt.Relaxation)
  % event-specific relaxation (full-space engine only)
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

% Choose the simulation engine
%===============================================================================
[Opt.Method,methodReason,seq] = sf_selectengine(Sys,Exp,Opt);
logmsg(1,'Simulation engine: %s %s',Opt.Method,methodReason);
subspaceEngine = strcmp(Opt.Method,'subspace');

if subspaceEngine
  Opt.fastSeq = seq;
  Opt.fastExp = Exp;
else
  % Full-space engine: set up events and indirect dimensions
  ExpProp = Exp;
  if isfield(ExpProp,'DetOperator')
    Opt.DetOperator = ExpProp.DetOperator;
    ExpProp = rmfield(ExpProp,'DetOperator');
  end
  if ischar(ExpProp.DetWindow) || isstring(ExpProp.DetWindow)
    if ~strcmp(ExpProp.DetWindow,'sequence')
      error('Exp.DetWindow must be a number, two numbers, or ''sequence''.');
    end
    ExpProp = rmfield(ExpProp,'DetWindow');
    ExpProp.DetSequence = 'all';
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
  [Events,Vary,Opt] = sf_buildevents(ExpProp,Opt);
  % Detected transients whose length changes along an indirect dimension
  % are not supported: sf_propagate stops at the first point where the
  % length changes
  if ~isempty(Vary)
    Vary.ConstantLength = true;
  end
  Opt.propEvents = Events;
  Opt.propVary = Vary;
  Opt.propPathways = sf_pathwayorders(Exp,Opt);
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
% Single-isotopologue simulation
%===============================================================================
% Exp and Opt have been parsed by saffron(), including the pulse sequence
% (Opt.fastSeq for the subspace engine, Opt.propEvents and Opt.propVary
% for the full-space engine).
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
if isfield(Sys,'nn') && any(Sys.nn(:)~=0) && subspaceEngine
  error('The subspace engine does not support nuclear-nuclear couplings (Sys.nn).');
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
% Electron coherence orders of the selected pathways (one row per
% pathway, NaN = any order), for the full-space engine. Empty: no
% selection.
function orders = sf_pathwayorders(Exp,Opt)
orders = [];
if ~isfield(Exp,'Pathways') || isempty(Exp.Pathways), return; end
pw = Exp.Pathways;
nPulses = sum(cellfun(@isstruct,Exp.Sequence));
if (ischar(pw) || isstring(pw)) && strcmp(pw,'all')
  return
end
if (ischar(pw) || isstring(pw)) && strcmp(pw,'echo')
  [seq,~] = sf_parsesequence(Exp,Opt);
  orders = unique(seq.co,'rows');
  return
end
if ischar(pw) || isstring(pw)
  pw = char(pw);
  if size(pw,2)~=nPulses
    error('Exp.Pathways has %d entries per pathway, but the sequence has %d pulses.',size(pw,2),nPulses);
  end
  if any(~ismember(pw(:),'+-0.'))
    error('Exp.Pathways can only contain ''+'', ''-'', ''0'' and ''.''.');
  end
  orders = nan(size(pw));
  orders(pw=='+') = +1;
  orders(pw=='-') = -1;
  orders(pw=='0') = 0;
else
  orders = pw;
  if size(orders,2)~=nPulses
    error('Exp.Pathways has %d entries per pathway, but the sequence has %d pulses.',size(orders,2),nPulses);
  end
end
end

%===============================================================================
% Dimensions whose first target is the RF carrier frequency
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
% Axes of the indirect dimensions: values of the first parameter listed in
% each Exp.DimN (its initial value plus the increments)
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
% Processing: baseline correction, apodization and Fourier transform of
% single-point signals along delay dimensions, giving info.f and info.fd;
% ENDOR spectra are returned as info.fd
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
% Plotting
%===============================================================================
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

% Labels of the nAx axes: indirect dimensions, then the transient time
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
