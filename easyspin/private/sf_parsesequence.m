% sf_parsesequence  Parse a user-defined experiment for the subspace engine
%
%   [sequence,reasons] = sf_parsesequence(Exp,idealPulses)
%   [sequence,reasons] = sf_parsesequence(Exp,idealPulses,relaxation)
%
%   Converts a user-defined experiment into the representation used by the
%   subspace engine (sf_engine_subspace). Each interval is a pulse followed
%   by a (possibly zero-length) delay.
%
%   Indirect dimensions are of three kinds:
%     1  linear increments of delays, computed by frequency-domain binning
%        (sf_peaks); dimensions are relabeled 1,2 in order of appearance
%     2  the RF carrier frequency (Exp.rfFreq) with a single RF pulse,
%        computed from the analytic response of each NMR transition
%     3  any other parameter of the pulses or delays, computed point by
%        point; also linear delay increments that do not fit a binning
%        scheme of sf_peaks (several incremented delays in one interval,
%        other combinations of incremented intervals)
%
%   Coherence transfer pathways are coded per interval as 1 = alpha
%   population, 2 = beta population, 3 = +1 coherence, 4 = -1 coherence.
%
%   Input:
%     Exp  experiment structure, with fields
%       Sequence         pulse sequence: cell array of pulses and delays (µs)
%       nPoints, Dim1, Dim2, ...   incrementation scheme
%       Pathways         coherence pathway selection, default 'all'
%       PhaseCycle       phase cycle of each pulse
%       DetPhase         detection phase (radians), default 0
%       rfFreq           RF carrier frequency (MHz), for RF pulses
%       mwFreq           microwave carrier frequency (GHz), optional
%       DetWindow, DetOperator, DetIntegrate, DetFreq, ResonatorFrequency,
%       ResonatorQL, FrequencyResponse, rfPhaseLocked
%                        only checked for unsupported features
%     idealPulses        true: treat finite pulses as ideal
%     relaxation         relaxation on/off per element of Exp.Sequence
%                        (Opt.Relaxation); a single value or omitted: on
%                        during all delays
%
%   Output:
%     seq  sequence structure (N intervals, nPathways pathways), with fields
%       N                number of intervals (number of pulses)
%       Pulse            1xN structure array, with fields Channel ('mw' or
%                        'rf'), Phase, tp, Frequency, Type, Par (the pulse
%                        as given), Flip, ideal, B1 (RF pulses)
%       t                delay after each pulse (µs), 1xN
%       tRelax           part of t with relaxation on (µs), 1xN
%       incRelax         true if the incremented delay of the interval
%                        has relaxation on, 1xN
%       delayInterval    interval of each delay in Exp.Sequence (0: before
%                        the first pulse)
%       isRF, nRF        true for RF pulses (1xN), number of RF pulses
%       rfFreq, rfAxis   RF carrier (MHz), and RF axis (MHz) of an RF sweep
%       ExciteWidth      excitation window (MHz) for idealPulses, or []
%       finitePulses     true if there are finite microwave pulses
%       DetPhase         detection phase (radians)
%       nPoints, kind    number of points and kind (1, 2, 3) of each dimension
%       mods             cell array with one element per dimension; kinds 1
%                        and 3: rows {'d' or 'p', delay or pulse index,
%                        pulse field, values}
%       Inc              relabeled kind-1 dimension of each interval (0:
%                        not incremented), 1xN
%       kind1Dims        original dimension of the relabeled dimensions
%       dt, nPointsInc   increments (µs) and points of the relabeled dimensions
%       IncSchemeID      incrementation scheme of sf_peaks
%       ctp              coherence transfer pathways, nPathways x N
%       co               electron coherence orders, nPathways x N
%       echoOrders       orders selected by Exp.Pathways = 'echo', before
%                        pathways with zero prefactor are removed (empty
%                        otherwise)
%       idxPulseL, idxPulseR   propagator blocks (1 = aa, 2 = bb, 3 = ab,
%                        4 = ba) of the pulses, nPathways x N
%       idxIncL, idxIncR manifolds (1 = alpha, 2 = beta) of the incremented
%                        free evolutions, one row per pathway
%       pathwayprefactor pulse and phase cycle prefactor of each pathway
%       nPathways        number of pathways
%       rfOn, rfOff      phase cycle weights with RF on and off
%     reasons   cell array of strings that list the features the subspace
%               engine does not support (empty if it supports all of them)

function [sequence,reasons] = sf_parsesequence(Exp,idealPulses,relaxation)

reasons = {};
sequence = struct;

Sequence = Exp.Sequence;
if ~iscell(Sequence)
  error('Exp.Sequence must be a cell array of pulses and delays.');
end

isPulse = cellfun(@isstruct,Sequence);
if ~any(isPulse)
  error('Exp.Sequence contains no pulses.');
end
for p = find(isPulse)
  if numel(Sequence{p})>1
    reasons{end+1} = 'simultaneous pulses';
  end
end
for p = find(~isPulse)
  if ~isnumeric(Sequence{p}) || ~isscalar(Sequence{p}) || Sequence{p}<0
    error('Delays in Exp.Sequence must be single nonnegative numbers (in µs).');
  end
end

% Intervals: pulse k and the delays after it. Delays before the first
% pulse have no effect on the initial state -Sz (no relaxation is applied
% before the first pulse). tRelax is the part of the delays of each
% interval with relaxation on (Opt.Relaxation per element of Exp.Sequence;
% entries for pulses have no effect, since there is no relaxation during
% pulses).
if nargin<3 || isempty(relaxation) || isscalar(relaxation)
  relaxed = true(1,numel(Sequence));
else
  relaxed = logical(relaxation(:).');
end
idxPulse = find(isPulse);
nPulses = numel(idxPulse);
idxDelay = find(~isPulse);
delayInterval = zeros(1,numel(idxDelay));  % interval of each delay (0: before first pulse)
t = zeros(1,nPulses);
tRelax = zeros(1,nPulses);
for iD = 1:numel(idxDelay)
  iInt = sum(idxPulse<idxDelay(iD));
  delayInterval(iD) = iInt;
  if iInt>0
    t(iInt) = t(iInt) + Sequence{idxDelay(iD)};
    if relaxed(idxDelay(iD))
      tRelax(iInt) = tRelax(iInt) + Sequence{idxDelay(iD)};
    end
  end
end
sequence.delayInterval = delayInterval;
sequence.t = t;
sequence.tRelax = tRelax;

% Pulses
Pulses = struct('Channel',cell(1,nPulses));
for k = 1:nPulses
  pulse_ = Sequence{idxPulse(k)};
  pulse_ = pulse_(1);
  if isfield(pulse_,'Channel') && ~isempty(pulse_.Channel)
    ch = lower(char(pulse_.Channel));
    if ~any(strcmp(ch,{'mw','rf'}))
      error('Pulse %d: Channel must be ''mw'' or ''rf''.',k);
    end
  else
    ch = 'mw';
  end
  Pulses(k).Channel = ch;
  if isfield(pulse_,'IQ')
    reasons{end+1} = 'user-defined pulse waveforms (IQ)'; %#ok<*AGROW>
  end
  if ~isfield(pulse_,'Phase') || isempty(pulse_.Phase), pulse_.Phase = 0; end
  Pulses(k).Phase = pulse_.Phase;
  if ~isfield(pulse_,'tp') || isempty(pulse_.tp), pulse_.tp = 0; end
  Pulses(k).tp = pulse_.tp;
  if ~isfield(pulse_,'Frequency') || isempty(pulse_.Frequency), pulse_.Frequency = 0; end
  Pulses(k).Frequency = pulse_.Frequency;
  if isfield(pulse_,'Type') && ~isempty(pulse_.Type)
    Pulses(k).Type = char(pulse_.Type);
  else
    Pulses(k).Type = 'rectangular';
  end
  Pulses(k).Par = pulse_;
  rect = any(strcmp(Pulses(k).Type,{'rectangular','rectangular/none'}));
  switch ch
    case 'mw'
      if ~rect
        reasons{end+1} = 'shaped or frequency-swept microwave pulses';
      end
      if numel(Pulses(k).Frequency)>1 && Pulses(k).Frequency(1)~=Pulses(k).Frequency(end)
        reasons{end+1} = 'frequency-swept microwave pulses';
      end
      Pulses(k).Frequency = Pulses(k).Frequency(1);
      if isfield(pulse_,'Flip') && ~isempty(pulse_.Flip)
        Pulses(k).Flip = pulse_.Flip;
      elseif isfield(pulse_,'Amplitude') && Pulses(k).tp>0
        Pulses(k).Flip = 2*pi*pulse_.Amplitude*Pulses(k).tp;
      elseif rect && ~isfield(pulse_,'IQ')
        error('Pulse %d: the flip angle (Flip) is missing.',k);
      else
        Pulses(k).Flip = pi;  % placeholder, shaped and IQ pulses are not supported
      end
      Pulses(k).ideal = Pulses(k).tp==0 || idealPulses;
    case 'rf'
      if Pulses(k).tp<=0
        error('Pulse %d: RF pulses need a pulse length (tp>0).',k);
      end
      if ~rect
        reasons{end+1} = 'shaped or frequency-swept RF pulses';
      end
      if numel(Pulses(k).Frequency)>1 && Pulses(k).Frequency(1)~=Pulses(k).Frequency(end)
        reasons{end+1} = 'frequency-swept RF pulses';
      end
      Pulses(k).Frequency = Pulses(k).Frequency(1);
      Pulses(k).ideal = false;
      Pulses(k).Flip = [];
      Pulses(k).B1 = [];
      if isfield(pulse_,'B1') && ~isempty(pulse_.B1)
        Pulses(k).B1 = pulse_.B1;
      elseif isfield(pulse_,'Flip') && ~isempty(pulse_.Flip)
        Pulses(k).Flip = pulse_.Flip;
      elseif ~rect
        Pulses(k).Flip = pi;  % placeholder, shaped RF pulses are not supported
      else
        error('Pulse %d: RF pulses need Flip (nominal) or B1 (in mT).',k);
      end
  end
end
sequence.Pulse = Pulses;
sequence.N = nPulses;

% RF pulses
isRF = strcmp({Pulses.Channel},'rf');
sequence.isRF = isRF;
sequence.nRFpulses = sum(isRF);
if sequence.nRFpulses>1
  reasons{end+1} = 'more than one RF pulse';
end
if sequence.nRFpulses>0 && ~isfield(Exp,'rfFreq')
  error('Exp.rfFreq (RF carrier frequency, in MHz) is required for RF pulses.');
end
if sequence.nRFpulses>0
  sequence.rfFreq = Exp.rfFreq;
end

% Excitation window for ideal-pulse orientation selection (idealPulses)
mwFinite = ~isRF & [Pulses.tp]>0;
if idealPulses && any(mwFinite)
  sequence.ExciteWidth = 1/min([Pulses(mwFinite).tp]);  % MHz
else
  sequence.ExciteWidth = [];
end
sequence.finitePulses = any(~[Pulses.ideal] & ~isRF);

% Without Exp.mwFreq, the subspace engine puts finite microwave pulses on
% resonance with every orientation (no orientation selection). A pulse
% frequency is then absolute (carrier at 0), which needs the full-space
% engine.
if ~isfield(Exp,'mwFreq')
  for p = find(~isRF & ~[Pulses.ideal])
    if any(Pulses(p).Par.Frequency~=0)
      reasons{end+1} = 'absolute pulse frequencies (no Exp.mwFreq)';
      break
    end
  end
end

% Detection
if isfield(Exp,'DetWindow') && ~(isnumeric(Exp.DetWindow) && all(Exp.DetWindow==0))
  reasons{end+1} = 'transient or shifted detection (Exp.DetWindow)';
end
if isfield(Exp,'DetOperator') && ~isempty(Exp.DetOperator)
  reasons{end+1} = 'detection operators (Exp.DetOperator)';
end
if isfield(Exp,'DetIntegrate') && Exp.DetIntegrate
  reasons{end+1} = 'echo integration (Exp.DetIntegrate)';
end
if isfield(Exp,'DetFreq') && ~(isfield(Exp,'mwFreq') && isequal(Exp.DetFreq,Exp.mwFreq))
  reasons{end+1} = 'Exp.DetFreq different from Exp.mwFreq';
end
if any(isfield(Exp,{'ResonatorFrequency','ResonatorQL','FrequencyResponse'}))
  reasons{end+1} = 'resonator';
end
if isfield(Exp,'DetPhase') && ~isempty(Exp.DetPhase)
  sequence.DetPhase = Exp.DetPhase;
else
  sequence.DetPhase = 0;
end

% Indirect dimensions
%-------------------------------------------------------------------------------
if isfield(Exp,'nPoints') && ~isempty(Exp.nPoints)
  nPoints = Exp.nPoints;
else
  nPoints = [];
end
nDims = numel(nPoints);
Inc = zeros(1,nPulses);       % kind-1 dimension label of each interval (before relabeling)
incRelax = true(1,nPulses);   % relaxation on for the incremented delay of each interval
dt = zeros(1,nDims);
kind = zeros(1,nDims);
mods = cell(1,nDims);   % kind 3: {target type, index, field, values}
dimTargets = cell(1,nDims);  % quantities changed by each dimension
dimValues = cell(1,nDims);   % their values per point
for d = 1:nDims
  dimField = sprintf('Dim%d',d);
  if ~isfield(Exp,dimField)
    error('Exp.%s is missing.',dimField);
  end
  Dim = Exp.(dimField);
  targets = {};
  values = {};
  for iLine = 1:size(Dim,1)
    tg = regexp(Dim{iLine,1},',','split');
    for q = 1:numel(tg)
      targets{end+1} = strtrim(tg{q});
      values{end+1} = Dim{iLine,2};
    end
  end
  % values per point for each target
  vals = cell(1,numel(targets));
  linStep = nan(1,numel(targets));
  for q = 1:numel(targets)
    v = values{q};
    if size(v,2)>1 && endsWith(targets{q},'.Frequency')
      % several values per point: a row of band-edge increments, or one
      % frequency band per row (as read by sf_buildevents); a column of
      % per-point offsets is handled below
      reasons{end+1} = 'indirect dimensions over frequency bands of pulses';
    elseif endsWith(targets{q},'.IQ')
      % one user-defined waveform per point, given in the pulse (no values);
      % the full-space reason is added below (pulse field 'IQ')
    elseif isscalar(v)
      vals{q} = v*(0:nPoints(d)-1);
      linStep(q) = v;
    else
      if numel(v)~=nPoints(d)
        error('Exp.%s: the number of values (%d) does not match Exp.nPoints(%d) = %d.',dimField,numel(v),d,nPoints(d));
      end
      vals{q} = v(:).';
      if nPoints(d)>1 && v(1)==0 && all(abs(diff(v,2))<=1e-12*max(abs(v)))
        linStep(q) = v(2);
      end
    end
  end
  isDelay = ~cellfun(@isempty,regexp(targets,'^d\d+$','once'));
  if all(isDelay) && linStep(1)>0 && all(linStep==linStep(1))
    % linear delay increments: binned (kind 1) if they fit a scheme of
    % sf_peaks (decided below), otherwise point by point (kind 3, with the
    % delay modifications in mods)
    kind(d) = 1;
    dt(d) = linStep(1);
    m = cell(numel(targets),4);
    for q = 1:numel(targets)
      iD = str2double(targets{q}(2:end));
      if iD<1 || iD>numel(idxDelay)
        error('Exp.%s: delay %s does not exist.',dimField,targets{q});
      end
      if delayInterval(iD)==0
        error('Exp.%s: incrementing a delay before the first pulse has no effect.',dimField);
      end
      m(q,:) = {'d',iD,'',vals{q}};
    end
    mods{d} = m;
  elseif all(strcmp(targets,'rfFreq'))
    kind(d) = 2;
    if isnan(linStep(1))
      sequence.rfAxis = Exp.rfFreq + vals{1};
    else
      sequence.rfAxis = Exp.rfFreq + linStep(1)*(0:nPoints(d)-1);
    end
  else
    kind(d) = 3;
    m = {};
    for q = 1:numel(targets)
      tg = targets{q};
      if isDelay(q)
        iD = str2double(tg(2:end));
        if iD<1 || iD>numel(idxDelay)
          error('Exp.%s: delay %s does not exist.',dimField,tg);
        end
        m(end+1,:) = {'d',iD,'',vals{q}};
      else
        tok = regexp(tg,'^p(\d+)\.(\w+)$','tokens','once');
        if isempty(tok)
          reasons{end+1} = sprintf('indirect dimension over ''%s''',tg);
          continue
        end
        iP = str2double(tok{1});
        if iP<1 || iP>nPulses
          error('Exp.%s: pulse p%d does not exist.',dimField,iP);
        end
        field = tok{2};
        if isRF(iP)
          reasons{end+1} = 'indirect dimensions over RF pulse parameters';
        elseif ~any(strcmp(field,{'Flip','Phase','tp','Frequency','Amplitude'}))
          reasons{end+1} = sprintf('indirect dimension over pulse field ''%s''',field);
        end
        m(end+1,:) = {'p',iP,field,vals{q}};
      end
    end
    mods{d} = m;
  end
  dimTargets{d} = targets;
  dimValues{d} = vals;
end
if sum(kind==2)>1 || (any(kind==2) && (any(kind==1) || any(kind==3)))
  reasons{end+1} = 'RF frequency sweeps combined with other indirect dimensions';
end
if any(kind==2) && sequence.nRFpulses~=1
  reasons{end+1} = 'RF frequency sweeps need exactly one RF pulse';
end
if sequence.nRFpulses>0 && ~any(kind==2)
  reasons{end+1} = 'RF pulses without an RF frequency sweep';
end
if isfield(Exp,'rfPhaseLocked') && Exp.rfPhaseLocked
  reasons{end+1} = 'phase-locked RF (Exp.rfPhaseLocked)';
end

% Binned delay increments: sf_peaks needs at most one incremented delay per
% interval, and the incremented intervals, labeled with their dimensions
% relabeled 1,2 in order of appearance, must form one of its schemes. The
% dimensions are added in the order of decreasing number of points, so that
% the larger ones stay binned. A dimension that does not fit is computed
% point by point (kind 3), which is exact but slower.
schemes = {[],1,[1 1],[1 2],[1 2 1],[1 1 2],[1 2 2 1],[1 1 2 2]};
schemeIDs = [0 1 2 11 12 14 15 17];
linDims = find(kind==1);
[~,bySize] = sort(nPoints(linDims),'descend');
for d = linDims(bySize)
  iD = [mods{d}{:,2}];
  intervals = delayInterval(iD);
  IncTrial = Inc;
  IncTrial(intervals) = d;
  [~,~,labels] = unique(IncTrial(IncTrial~=0),'stable');
  fits = numel(unique(intervals))==numel(intervals) && all(Inc(intervals)==0) && ...
    any(cellfun(@(s)isequal(s,labels.'),schemes));
  if fits
    Inc = IncTrial;
    incRelax(intervals) = relaxed(idxDelay(iD));
  else
    kind(d) = 3;
  end
end
sequence.nPoints = nPoints;
sequence.kind = kind;
sequence.mods = mods;

% Relabel kind-1 dimensions 1,2 in order of appearance, determine the
% incrementation scheme ID used by sf_peaks
incDims = Inc(Inc~=0);
order = unique(incDims,'stable');
relabel = zeros(1,nDims);
relabel(order) = 1:numel(order);
IncScheme = relabel(incDims);
sequence.Inc = zeros(1,nPulses);
sequence.Inc(Inc~=0) = IncScheme;
sequence.incRelax = incRelax;
sequence.kind1Dims = order;              % original dimension of relabeled dims 1,2
sequence.dt = dt(order);
sequence.nPointsInc = nPoints(order);
if isempty(IncScheme)
  sequence.IncSchemeID = 0;
else
  sequence.IncSchemeID = schemeIDs(cellfun(@(s)isequal(s,IncScheme),schemes));
end

% Parse/generate coherence transfer pathways
if isfield(Exp,'Pathways') && ~isempty(Exp.Pathways)
  pathways = Exp.Pathways;
else
  pathways = 'all';
end
[ctp,co,sequence.echoOrders] = selectPathways(pathways,sequence,dimTargets,dimValues);
nPathways = size(ctp,1);

% Block information for left and right propagation
[idxFreeL,idxFreeR,idxPulseL,idxPulseR] = pathwayparser(ctp);

% Pulse transfer prefactors (ideal mw pulses) and phase cycle weights
prefactor = ones(nPathways,1);
rfOn = 1; rfOff = 0;
PC = {};
if isfield(Exp,'PhaseCycle') && ~isempty(Exp.PhaseCycle)
  PC = Exp.PhaseCycle;
  if ~iscell(PC), PC = {PC}; end
  if numel(PC)>nPulses
    error('Exp.PhaseCycle has more entries (%d) than there are pulses (%d).',numel(PC),nPulses);
  end
end
PC(end+1:nPulses) = {[]};
PC(cellfun(@isempty,PC)) = {[0 1 1]};
dco = diff([zeros(nPathways,1) co],1,2);  % change in coherence order at each pulse
for p = 1:nPulses
  PCstep = PC{p};
  if size(PCstep,2)==2, PCstep(:,3) = 1; end
  if size(PCstep,2)~=3
    error('Exp.PhaseCycle{%d} must have 2 columns [phase weight] or 3 columns [phase weight scale].',p);
  end
  cycPhase = PCstep(:,1);
  cycWeight = PCstep(:,2)/sum(abs(PCstep(:,2)));
  cycOnOff = PCstep(:,3);

  if isRF(p)
    on = cycOnOff~=0;
    if any(cycOnOff(on)~=1)
      reasons{end+1} = 'RF amplitude scales other than 0 and 1 in phase cycles';
    end
    rfOn = sum(cycWeight(on));
    rfOff = sum(cycWeight(~on));
  elseif Pulses(p).ideal
    theta = Pulses(p).Flip*cycOnOff;
    phi = Pulses(p).Phase + cycPhase;
    c = cos(theta/2);
    sn = sin(theta/2);
    U = [c c -1i*sn.*exp(-1i*phi) -1i*sn.*exp(1i*phi)];  % blocks aa, bb, ab, ba; one row per cycle step
    f = (U(:,idxPulseL(:,p)).*conj(U(:,idxPulseR(:,p)))).'*cycWeight;
    prefactor = prefactor.*f;
  else
    if any(cycOnOff~=1)
      reasons{end+1} = 'amplitude scales in phase cycles of finite microwave pulses';
    end
    f = exp(-1i*dco(:,p)*cycPhase.')*cycWeight;
    prefactor = prefactor.*f;
  end
end

% Remove pathways with negligible prefactor
rmv = abs(prefactor)<1e-10;
if all(rmv)
  % no signal (e.g. a pi pulse creates no coherence): keep one pathway
  % with zero amplitude
  rmv(1) = false;
  prefactor(1) = 0;
end
if any(rmv)
  ctp(rmv,:) = [];
  co(rmv,:) = [];
  idxFreeL(rmv,:) = [];
  idxFreeR(rmv,:) = [];
  idxPulseL(rmv,:) = [];
  idxPulseR(rmv,:) = [];
  prefactor(rmv) = [];
end

sequence.ctp = ctp;
sequence.co = co;
sequence.idxPulseL = idxPulseL;
sequence.idxPulseR = idxPulseR;
sequence.idxIncL = idxFreeL(:,sequence.Inc~=0);
sequence.idxIncR = idxFreeR(:,sequence.Inc~=0);
% order the incremented-interval free indices by sequence order (as sf_peaks)
sequence.pathwayprefactor = prefactor;
sequence.nPathways = numel(prefactor);
sequence.rfOn = rfOn;
sequence.rfOff = rfOff;

reasons = unique(reasons,'stable');

end

%-------------------------------------------------------------------------------
% selectPathways  Coherence transfer pathways selected by Exp.Pathways
%
%   [ctp,co,echoOrders] = selectPathways(pathways,sequence,dimTargets,dimValues)
%
%   Enumerates the pathways that can be detected and applies Exp.Pathways.
%   With 'echo', selects the pathways that refocus at the detection point at
%   all points of the experiment: at least one coherence order inversion,
%   zero net coherence evolution at the first point, and zero change of it
%   at every point along each dimension. The contributions of the dimensions
%   add up, so this holds for the entire grid.
%
%   Input:
%     pathways    pathway selection (Exp.Pathways): 'all', 'echo', a character
%                 array of orders per pulse ('+', '-', '0', 'c' = +1 or -1,
%                 '.' = any), or a numeric matrix (NaN = any)
%     sequence    parsed sequence; uses N, t, isRF, Pulse, delayInterval, nPoints
%     dimTargets  quantities changed by each dimension, cell array with one
%                 cell array of strings per dimension
%     dimValues   values per point of these quantities, in the same layout
%
%   Output:
%     ctp         coherence transfer pathways, nPathways x nPulses
%     co          electron coherence orders, nPathways x nPulses
%     echoOrders  co for Exp.Pathways = 'echo', otherwise empty

function [ctp,co,echoOrders] = selectPathways(pathways,sequence,dimTargets,dimValues)

nPulses = sequence.N;
isRF = sequence.isRF;  % whether a pulse is an RF pulse

% Enumerate all pathways that end in electron coherence order -1 (S-, detected
% by S+). Encoding: 1=aa (Salpha), 2=bb (Sbeta), 3=ab (S+), 4=ba (S-).
ctp = 2;  % Sbeta (= identity/2-Sz) before first pulse
for p = 1:nPulses-1  % loop over pulses
  ctp = [repelem(ctp,4,1) repmat((1:4).',size(ctp,1),1)];  % each pulse gives 4 new states
end
ctp(:,nPulses+1) = 4;  % S- after last pulse

% RF pulses don't change the electron state
keep = all(ctp(:,[false isRF])==ctp(:,[isRF false]),2);
ctp = ctp(keep,2:end);  % also remove initial state

% Determine coherence orders for all intervals
co = (ctp==3) - (ctp==4);  % +1 and -1 coherence, 0 for populations

echoOrders = [];

if strcmp(pathways,'all')
  selected = true(size(co,1),1);

elseif strcmp(pathways,'echo')
  % Evolution time of each interval at the first point, and its change at
  % each point along each dimension. The spins evolve during RF pulses, so their
  % lengths count as evolution time (the lengths of microwave pulses do not).
  tFirst = sequence.t + isRF.*[sequence.Pulse.tp];
  changes = zeros(nPulses,0);
  for d = 1:numel(dimTargets)
    change = zeros(nPulses,sequence.nPoints(d));
    for k = 1:numel(dimTargets{d})
      target = dimTargets{d}{k};
      pulseField = regexp(target,'^p(\d+)\.(tp?|Position)$','tokens','once');
      signs = +1;
      if ~isempty(regexp(target,'^d\d+$','once'))
        intervals = sequence.delayInterval(str2double(target(2:end)));
      elseif ~isempty(pulseField) && strcmp(pulseField{2},'Position')
        % a moved pulse lengthens the delay before it and shortens the one after it
        intervals = str2double(pulseField{1}) - [1 0];
        signs = [+1 -1];
      elseif ~isempty(pulseField) && isRF(str2double(pulseField{1}))
        intervals = str2double(pulseField{1});
      else
        continue
      end
      signs = signs(intervals>0);  % delays before the first pulse have no effect
      intervals = intervals(intervals>0);
      if isempty(intervals), continue; end
      v = dimValues{d}{k};
      tFirst(intervals) = tFirst(intervals) + signs*v(1);
      change(intervals,:) = change(intervals,:) + signs.'*(v-v(1));
    end
    changes = [changes change];
  end

  % Determine echo pathways with echoes that refocus at the end of the sequence
  echo = any(co==+1,2);
  times = [tFirst.' changes];
  threshold = 1e-12*max([1; abs(times(:))]);
  selected = echo & all(abs(co*times)<threshold,2);
  if ~any(selected)
    error('No pathway refocuses at the detection point at all points of this experiment.');
  end
  echoOrders = co(selected,:);

else
  if size(pathways,2)~=nPulses
    error('Exp.Pathways has %d entries per pathway, but the sequence has %d pulses.',size(pathways,2),nPulses);
  end
  orders = pathways;
  if ischar(pathways)
    if any(~ismember(pathways(:),'+-0c.'))
      error('Exp.Pathways can only contain ''+'', ''-'', ''0'', ''c'' and ''.''.');
    end
    orders = nan(size(pathways));  % NaN: any order
    orders(pathways=='+') = +1;
    orders(pathways=='-') = -1;
    orders(pathways=='0') = 0;
    orders(pathways=='c') = Inf;  % Inf: +1 or -1
  end
  before = [zeros(size(orders,1),1) orders(:,1:end-1)];
  change = abs(orders-before);  % NaN or Inf if an order is not fixed
  k = find(isRF & any(change>0 & change<Inf,1),1);
  if ~isempty(k)
    error('Exp.Pathways: RF pulses (pulse %d) cannot change the electron coherence order.',k);
  end
  selected = false(size(co,1),1);
  for r = 1:size(orders,1)
    o = orders(r,:);
    selected = selected | all(co==o | isnan(o) | (isinf(o) & co~=0),2);
  end
end

ctp = ctp(selected,:);
co = co(selected,:);
if isempty(ctp)
  error('No coherence transfer pathway is detectable (ending in -1) with Exp.Pathways.');
end

end

%-------------------------------------------------------------------------------
% pathwayparser  Propagator blocks of coherence transfer pathways
%
%   [freeL,freeR,pulseL,pulseR] = pathwayparser(ctp)
%
%   Input:
%     ctp             coherence transfer pathways, nPathways x nPulses array, with
%                     1 = Salpha, 2 = Sbeta, 3 = S+, 4 = S-
%
%   Output:
%     freeL, freeR    left and right manifolds (1 = alpha, 2 = beta) of each
%                     free evolution, nPathways x nPulses array
%     pulseL, pulseR  left and right propagator blocks (1 = aa, 2 = bb, 3 =
%                     ab, 4 = ba) of each pulse, nPathways x nPulses array

function [freeL,freeR,pulseL,pulseR] = pathwayparser(ctp)

% General superpropagator
%      Sa  Sb  S+  S-
% Sa   aa  ++  a+  +a
% Sb   --  bb  -b  b-
% S+   a-  +b  ab  +-
% S-   -a  b+  -+  ba
propL = [1 3 1 3; 4 2 4 2; 1 3 1 3; 4 2 4 2];
propR = [1 3 3 1; 4 2 2 4; 4 2 2 4; 1 3 3 1];

before = [ones(size(ctp,1),1)*2 ctp(:,1:end-1)];
idx = ctp + (before-1)*4;
pulseL = propL(idx);
pulseR = propR(idx);

idx = 5*ctp - 4;
freeL = propL(idx);
freeR = propR(idx);

end
