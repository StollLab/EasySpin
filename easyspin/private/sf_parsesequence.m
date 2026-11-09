% sf_parsesequence  Parse a user-defined experiment for the subspace engine
%
%   [seq,reasons] = sf_parsesequence(Exp,Opt)
%
%   Converts a user-defined experiment into the representation used by the
%   subspace engine (sf_engine_subspace). Each interval is a pulse followed
%   by a (possibly zero) delay.
%
%   Indirect dimensions are of three kinds:
%     1  linear increments of delays, computed by frequency-domain binning
%        (sf_peaks); dimensions are relabeled 1,2 in order of appearance
%     2  the RF carrier frequency (Exp.rfFreq) with a single RF pulse,
%        computed from the analytic response of each NMR transition
%     3  any other parameter of the pulses or delays, computed point by point
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
%     Opt  options structure, with fields
%       IdealPulses      true: treat finite pulses as ideal
%
%   Output:
%     seq  sequence structure (N intervals, nPathways pathways), with fields
%       N                number of intervals (number of pulses)
%       Pulse            1xN structure array, with fields Channel ('mw' or
%                        'rf'), Phase, tp, Frequency, Type, Par (the pulse
%                        as given), Flip, ideal, B1 (RF pulses)
%       t                delay after each pulse (µs), 1xN
%       delayInterval    interval of each delay in Exp.Sequence (0: before
%                        the first pulse)
%       isRF, nRF        true for RF pulses (1xN), number of RF pulses
%       rfFreq, rfAxis   RF carrier (MHz), and RF axis (MHz) of an RF sweep
%       ExciteWidth      excitation window (MHz) for Opt.IdealPulses, or []
%       finitePulses     true if there are finite microwave pulses
%       DetPhase         detection phase (radians)
%       nPoints, kind    number of points and kind (1, 2, 3) of each dimension
%       mods             cell array with one element per dimension; kind 3:
%                        rows {'d' or 'p', delay or pulse index, pulse
%                        field, values}
%       Inc              relabeled kind-1 dimension of each interval (0:
%                        not incremented), 1xN
%       kind1Dims        original dimension of the relabeled dimensions
%       dt, nPointsInc   increments (µs) and points of the relabeled dimensions
%       IncSchemeID      incrementation scheme of sf_peaks, NaN if unsupported
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

function [seq,reasons] = sf_parsesequence(Exp,Opt)

reasons = {};
seq = struct;

Sequence = Exp.Sequence;
if ~iscell(Sequence)
  error('Exp.Sequence must be a cell array of pulses and delays.');
end

isPulse = cellfun(@isstruct,Sequence);
if ~any(isPulse)
  error('Exp.Sequence contains no pulses.');
end
for k = find(isPulse)
  if numel(Sequence{k})>1
    reasons{end+1} = 'simultaneous pulses';
  end
end
for k = find(~isPulse)
  if ~isnumeric(Sequence{k}) || ~isscalar(Sequence{k}) || Sequence{k}<0
    error('Delays in Exp.Sequence must be single nonnegative numbers (in µs).');
  end
end

% Intervals: pulse k and the delays after it. Delays before the first
% pulse have no effect on the initial state -Sz (no relaxation is applied
% before the first pulse).
pulseIdx = find(isPulse);
N = numel(pulseIdx);
delayIdx = find(~isPulse);
delayInterval = zeros(1,numel(delayIdx));  % interval of each delay (0: before first pulse)
t = zeros(1,N);
for iD = 1:numel(delayIdx)
  iInt = sum(pulseIdx<delayIdx(iD));
  delayInterval(iD) = iInt;
  if iInt>0
    t(iInt) = t(iInt) + Sequence{delayIdx(iD)};
  end
end
seq.delayInterval = delayInterval;

% Pulses
idealPulses = isfield(Opt,'IdealPulses') && Opt.IdealPulses;
P = struct('Channel',cell(1,N));
for k = 1:N
  p_ = Sequence{pulseIdx(k)};
  p_ = p_(1);
  if isfield(p_,'Channel') && ~isempty(p_.Channel)
    ch = lower(char(p_.Channel));
    if ~any(strcmp(ch,{'mw','rf'}))
      error('Pulse %d: Channel must be ''mw'' or ''rf''.',k);
    end
  else
    ch = 'mw';
  end
  P(k).Channel = ch;
  if isfield(p_,'IQ')
    reasons{end+1} = 'user-defined pulse waveforms (IQ)'; %#ok<*AGROW>
  end
  if ~isfield(p_,'Phase') || isempty(p_.Phase), p_.Phase = 0; end
  P(k).Phase = p_.Phase;
  if ~isfield(p_,'tp') || isempty(p_.tp), p_.tp = 0; end
  P(k).tp = p_.tp;
  if ~isfield(p_,'Frequency') || isempty(p_.Frequency), p_.Frequency = 0; end
  P(k).Frequency = p_.Frequency;
  if isfield(p_,'Type') && ~isempty(p_.Type)
    P(k).Type = char(p_.Type);
  else
    P(k).Type = 'rectangular';
  end
  P(k).Par = p_;
  rect = any(strcmp(P(k).Type,{'rectangular','rectangular/none'}));
  switch ch
    case 'mw'
      if ~rect
        reasons{end+1} = 'shaped or frequency-swept microwave pulses';
      end
      if numel(P(k).Frequency)>1 && P(k).Frequency(1)~=P(k).Frequency(end)
        reasons{end+1} = 'frequency-swept microwave pulses';
      end
      P(k).Frequency = P(k).Frequency(1);
      if isfield(p_,'Flip') && ~isempty(p_.Flip)
        P(k).Flip = p_.Flip;
      elseif isfield(p_,'Amplitude') && P(k).tp>0
        P(k).Flip = 2*pi*p_.Amplitude*P(k).tp;
      elseif rect && ~isfield(p_,'IQ')
        error('Pulse %d: the flip angle (Flip) is missing.',k);
      else
        P(k).Flip = pi;  % placeholder, shaped and IQ pulses are not supported
      end
      P(k).ideal = P(k).tp==0 || idealPulses;
    case 'rf'
      if P(k).tp<=0
        error('Pulse %d: RF pulses need a pulse length (tp>0).',k);
      end
      if ~rect
        reasons{end+1} = 'shaped or frequency-swept RF pulses';
      end
      if numel(P(k).Frequency)>1 && P(k).Frequency(1)~=P(k).Frequency(end)
        reasons{end+1} = 'frequency-swept RF pulses';
      end
      P(k).Frequency = P(k).Frequency(1);
      P(k).ideal = false;
      P(k).Flip = [];
      P(k).B1 = [];
      if isfield(p_,'B1') && ~isempty(p_.B1)
        P(k).B1 = p_.B1;
      elseif isfield(p_,'Flip') && ~isempty(p_.Flip)
        P(k).Flip = p_.Flip;
      elseif ~rect
        P(k).Flip = pi;  % placeholder, shaped RF pulses are not supported
      else
        error('Pulse %d: RF pulses need Flip (nominal) or B1 (in mT).',k);
      end
  end
end
seq.Pulse = P;
seq.N = N;
seq.t = t;
isRF = strcmp({P.Channel},'rf');
seq.isRF = isRF;
seq.nRF = sum(isRF);
if seq.nRF>1
  reasons{end+1} = 'more than one RF pulse';
end
if seq.nRF>0 && ~isfield(Exp,'rfFreq')
  error('Exp.rfFreq (RF carrier frequency, in MHz) is required for RF pulses.');
end
if seq.nRF>0
  seq.rfFreq = Exp.rfFreq;
end

% Excitation window for ideal-pulse orientation selection (Opt.IdealPulses)
mwFinite = ~isRF & [P.tp]>0;
if idealPulses && any(mwFinite)
  seq.ExciteWidth = 1/min([P(mwFinite).tp]);  % MHz
else
  seq.ExciteWidth = [];
end
seq.finitePulses = any(~[P.ideal] & ~isRF);

% Without Exp.mwFreq, the subspace engine puts finite microwave pulses on
% resonance with every orientation (no orientation selection). A pulse
% frequency is then absolute (carrier at 0), which needs the full-space
% engine.
if ~isfield(Exp,'mwFreq')
  for k = find(~isRF & ~[P.ideal])
    if any(P(k).Par.Frequency~=0)
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
  seq.DetPhase = Exp.DetPhase;
else
  seq.DetPhase = 0;
end

% Indirect dimensions
%-------------------------------------------------------------------------------
if isfield(Exp,'nPoints') && ~isempty(Exp.nPoints)
  nPoints = Exp.nPoints;
else
  nPoints = [];
end
nDims = numel(nPoints);
Inc = zeros(1,N);       % kind-1 dimension label of each interval (before relabeling)
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
    for k = 1:numel(tg)
      targets{end+1} = strtrim(tg{k});
      values{end+1} = Dim{iLine,2};
    end
  end
  % values per point for each target
  vals = cell(1,numel(targets));
  linStep = nan(1,numel(targets));
  for k = 1:numel(targets)
    v = values{k};
    if size(v,2)>1 && ~isempty(regexp(targets{k},'\.Frequency$','once'))
      % several values per point: a row of band-edge increments, or one
      % frequency band per row (as read by sf_buildevents); a column of
      % per-point offsets is handled below
      reasons{end+1} = 'indirect dimensions over frequency bands of pulses';
      vals{k} = [];
      continue
    end
    if ~isempty(regexp(targets{k},'\.IQ$','once'))
      % one user-defined waveform per point, given in the pulse (no values);
      % the full-space reason is added below (pulse field 'IQ')
      vals{k} = [];
      continue
    end
    if isscalar(v)
      vals{k} = v*(0:nPoints(d)-1);
      linStep(k) = v;
    else
      if numel(v)~=nPoints(d)
        error('Exp.%s: the number of values (%d) does not match Exp.nPoints(%d) = %d.',dimField,numel(v),d,nPoints(d));
      end
      vals{k} = v(:).';
      if nPoints(d)>1 && v(1)==0 && all(abs(diff(v,2))<=1e-12*max(abs(v)))
        linStep(k) = v(2);
      end
    end
  end
  isDelay = ~cellfun(@isempty,regexp(targets,'^d\d+$','once'));
  isRFfreq = strcmp(targets,'rfFreq');
  if all(isDelay) && all(~isnan(linStep)) && numel(unique(linStep))==1 && linStep(1)>0
    kind(d) = 1;
    dt(d) = linStep(1);
    for k = 1:numel(targets)
      iD = str2double(targets{k}(2:end));
      if iD<1 || iD>numel(delayIdx)
        error('Exp.%s: delay %s does not exist.',dimField,targets{k});
      end
      iInt = delayInterval(iD);
      if iInt==0
        error('Exp.%s: incrementing a delay before the first pulse has no effect.',dimField);
      end
      if Inc(iInt)~=0
        reasons{end+1} = 'several incremented delays in one interval';
      end
      Inc(iInt) = d;
    end
  elseif all(isRFfreq)
    kind(d) = 2;
    if isnan(linStep(1))
      seq.rfAxis = Exp.rfFreq + vals{1};
    else
      seq.rfAxis = Exp.rfFreq + linStep(1)*(0:nPoints(d)-1);
    end
  else
    kind(d) = 3;
    m = {};
    for k = 1:numel(targets)
      tg = targets{k};
      if isDelay(k)
        iD = str2double(tg(2:end));
        if iD<1 || iD>numel(delayIdx)
          error('Exp.%s: delay %s does not exist.',dimField,tg);
        end
        m(end+1,:) = {'d',iD,'',vals{k}};
      else
        tok = regexp(tg,'^p(\d+)\.(\w+)$','tokens','once');
        if isempty(tok)
          reasons{end+1} = sprintf('indirect dimension over ''%s''',tg);
          continue
        end
        iP = str2double(tok{1});
        if iP<1 || iP>N
          error('Exp.%s: pulse p%d does not exist.',dimField,iP);
        end
        field = tok{2};
        if isRF(iP)
          reasons{end+1} = 'indirect dimensions over RF pulse parameters';
        elseif ~any(strcmp(field,{'Flip','Phase','tp','Frequency','Amplitude'}))
          reasons{end+1} = sprintf('indirect dimension over pulse field ''%s''',field);
        end
        m(end+1,:) = {'p',iP,field,vals{k}};
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
if any(kind==2) && seq.nRF~=1
  reasons{end+1} = 'RF frequency sweeps need exactly one RF pulse';
end
if seq.nRF>0 && ~any(kind==2)
  reasons{end+1} = 'RF pulses without an RF frequency sweep';
end
if isfield(Exp,'rfPhaseLocked') && Exp.rfPhaseLocked
  reasons{end+1} = 'phase-locked RF (Exp.rfPhaseLocked)';
end
seq.nPoints = nPoints;
seq.kind = kind;
seq.mods = mods;

% Relabel kind-1 dimensions 1,2 in order of appearance, determine the
% incrementation scheme ID used by sf_peaks
incDims = Inc(Inc~=0);
order = unique(incDims,'stable');
relabel = zeros(1,nDims);
relabel(order) = 1:numel(order);
IncScheme = relabel(incDims);
seq.Inc = zeros(1,N);
seq.Inc(Inc~=0) = IncScheme;
seq.kind1Dims = order;              % original dimension of relabeled dims 1,2
seq.dt = dt(order);
seq.nPointsInc = nPoints(order);
schemes = {[],1,[1 1],[1 2],[1 2 1],[1 1 2],[1 2 2 1],[1 1 2 2]};
schemeIDs = [0 1 2 11 12 14 15 17];
if isempty(IncScheme)
  idx = 1;
else
  idx = find(cellfun(@(s)isequal(s,IncScheme),schemes));
end
if isempty(idx)
  reasons{end+1} = 'this combination of incremented delays';
  seq.IncSchemeID = NaN;
else
  seq.IncSchemeID = schemeIDs(idx);
end

% Coherence transfer pathways
if isfield(Exp,'Pathways') && ~isempty(Exp.Pathways)
  pathways = Exp.Pathways;
else
  pathways = 'all';
end
if isstring(pathways), pathways = char(pathways); end
[ctp,co,seq.echoOrders] = selectPathways(pathways,seq,dimTargets,dimValues);

% Pulse transfer prefactors (ideal mw pulses) and phase cycle weights
[idxFreeL,idxFreeR,idxPulseL,idxPulseR] = pathwayparser(ctp);
nPathways = size(ctp,1);
pref = ones(nPathways,1);
rfOn = 1; rfOff = 0;
PC = {};
if isfield(Exp,'PhaseCycle') && ~isempty(Exp.PhaseCycle)
  PC = Exp.PhaseCycle;
  if ~iscell(PC), PC = {PC}; end
  if numel(PC)>N
    error('Exp.PhaseCycle has more entries (%d) than there are pulses (%d).',numel(PC),N);
  end
end
coPrev = [zeros(nPathways,1) co(:,1:end-1)];
for k = 1:N
  if k<=numel(PC) && ~isempty(PC{k})
    cyc = PC{k};
    if size(cyc,2)==2, cyc(:,3) = 1; end
    if size(cyc,2)~=3
      error('Exp.PhaseCycle{%d} must have 2 columns [phase weight] or 3 columns [phase weight scale].',k);
    end
  else
    cyc = [0 1 1];
  end
  norm_ = sum(abs(cyc(:,2)));
  if isRF(k)
    on = cyc(:,3)~=0;
    if any(cyc(on,3)~=1)
      reasons{end+1} = 'RF amplitude scales other than 0 and 1 in phase cycles';
    end
    rfOn = sum(cyc(on,2))/norm_;
    rfOff = sum(cyc(~on,2))/norm_;
    continue
  end
  if P(k).ideal
    f = zeros(nPathways,1);
    for s = 1:size(cyc,1)
      theta = P(k).Flip*cyc(s,3);
      phi = P(k).Phase + cyc(s,1);
      c = cos(theta/2);
      sn = sin(theta/2);
      pL = [c c -1i*sn*exp(-1i*phi) -1i*sn*exp(+1i*phi)];
      pR = [c c +1i*sn*exp(+1i*phi) +1i*sn*exp(-1i*phi)];
      f = f + cyc(s,2)*pL(idxPulseL(:,k)).'.*pR(idxPulseR(:,k)).';
    end
    pref = pref.*f/norm_;
  else
    if any(cyc(:,3)~=1)
      reasons{end+1} = 'amplitude scales in phase cycles of finite microwave pulses';
    end
    dp = co(:,k)-coPrev(:,k);
    f = exp(-1i*dp*cyc(:,1).')*cyc(:,2);
    pref = pref.*f/norm_;
  end
end
rmv = abs(pref)<1e-10;
if all(rmv)
  % no signal (e.g. a pi pulse creates no coherence): keep one pathway
  % with zero amplitude
  rmv(1) = false;
  pref(1) = 0;
end
ctp(rmv,:) = []; co(rmv,:) = [];
idxFreeL(rmv,:) = []; idxFreeR(rmv,:) = [];
idxPulseL(rmv,:) = []; idxPulseR(rmv,:) = [];
pref(rmv) = [];

seq.ctp = ctp;
seq.co = co;
seq.idxPulseL = idxPulseL;
seq.idxPulseR = idxPulseR;
seq.idxIncL = idxFreeL(:,seq.Inc~=0);
seq.idxIncR = idxFreeR(:,seq.Inc~=0);
% order the incremented-interval free indices by sequence order (as sf_peaks)
seq.pathwayprefactor = pref;
seq.nPathways = numel(pref);
seq.rfOn = rfOn;
seq.rfOff = rfOff;

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
%   [FreeL,FreeR,PulsL,PulsR] = pathwayparser(ctp)
%
%   Input:
%     ctp           coherence transfer pathways, nPathways x N, with 1 =
%                   alpha population, 2 = beta population, 3 = +1 coherence,
%                   4 = -1 coherence
%
%   Output:
%     FreeL, FreeR  left and right manifolds (1 = alpha, 2 = beta) of each
%                   free evolution, nPathways x N
%     PulsL, PulsR  left and right propagator blocks (1 = aa, 2 = bb, 3 =
%                   ab, 4 = ba) of each pulse, nPathways x N
function [FreeL,FreeR,PulsL,PulsR] = pathwayparser(ctp)

% General superpropagator
%      Sa  Sb  S+  S-
% Sa   aa  ++  a+  +a
% Sb   --  bb  -b  b-
% S+   a-  +b  ab  +-
% S-   -a  b+  -+  ba
PropL = [1 3 1 3; 4 2 4 2; 1 3 1 3; 4 2 4 2];
PropR = [1 3 3 1; 4 2 2 4; 4 2 2 4; 1 3 3 1];

before = [ones(size(ctp,1),1)*2 ctp(:,1:end-1)];
idx = ctp + (before-1)*4;
PulsL = PropL(idx);
PulsR = PropR(idx);

idx = 5*ctp - 4;
FreeL = PropL(idx);
FreeR = PropR(idx);

end
