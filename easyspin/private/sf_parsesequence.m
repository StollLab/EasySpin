% sf_parsesequence  Parse a user-defined experiment for the subspace engine
%
%   [seq,reasons] = sf_parsesequence(Exp,Opt)
%
% Converts a user-defined experiment (the pulse sequence in Exp.Sequence as
% a cell array of pulses and delays, the incrementation scheme in
% Exp.nPoints, Exp.Dim1, ..., and Exp.Pathways, Exp.PhaseCycle,
% Exp.DetPhase) into the representation used by the subspace engine
% (sf_engine_subspace). Each interval is a pulse followed by a (possibly zero)
% delay. reasons is a cell array of strings that list the features the
% subspace engine does not support (empty if it supports all of them).
%
% Indirect dimensions are of three kinds:
%   1  linear increments of delays, computed by frequency-domain binning
%      (sf_peaks); dimensions are relabeled 1,2 in order of appearance
%   2  the RF carrier frequency (Exp.rfFreq) with a single RF pulse,
%      computed from the analytic response of each NMR transition
%   3  any other parameter of the pulses or delays, computed point by point
%
% Coherence transfer pathways are coded per interval as 1 = alpha
% population, 2 = beta population, 3 = +1 coherence, 4 = -1 coherence.

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
%-------------------------------------------------------------------------------
% All pathways that end in electron coherence order -1 (detected by S+).
% RF pulses don't change the electron state. The initial state is the
% beta population (the -Sz state differs from it by the identity/2, which
% gives no signal).
coherenceOrder = [0 0 +1 -1];
if N>1
  [idx_{1:N-1}] = ndgrid(1:4);
  ctp = fliplr(reshape(cat(N,idx_{:}),[],N-1));
  ctp(:,N) = 4;
else
  ctp = 4;
end
keep = true(size(ctp,1),1);
for k = find(isRF)
  if k==1
    before = 2*ones(size(ctp,1),1);
  else
    before = ctp(:,k-1);
  end
  keep = keep & ctp(:,k)==before;
end
ctp = ctp(keep,:);
co = coherenceOrder(ctp);

if ~isfield(Exp,'Pathways') || isempty(Exp.Pathways)
  Exp.Pathways = 'all';
end
pw = Exp.Pathways;
if (ischar(pw) || isstring(pw)) && any(strcmp(pw,{'all','echo'}))
  if strcmp(pw,'echo')
    % the spins evolve during RF pulses, so their lengths count as
    % evolution time (the lengths of microwave pulses do not)
    tEvol = t;
    tEvol(isRF) = tEvol(isRF) + [P(isRF).tp];
    sel = echoPathways(co,tEvol,seq.Inc);
    ctp = ctp(sel,:);
    co = co(sel,:);
    if isempty(ctp)
      error('No pathway refocuses at the detection point with this sequence and these delays.');
    end
  end
else
  if ischar(pw) || isstring(pw)
    pw = char(pw);
    if size(pw,2)~=N
      error('Exp.Pathways has %d entries per pathway, but the sequence has %d pulses.',size(pw,2),N);
    end
    sel_ = nan(size(pw));
    sel_(pw=='+') = +1;
    sel_(pw=='-') = -1;
    sel_(pw=='0') = 0;
    bad = ~ismember(pw,'+-0.');
    if any(bad(:))
      error('Exp.Pathways can only contain ''+'', ''-'', ''0'' and ''.''.');
    end
  else
    sel_ = pw;
    if size(sel_,2)~=N
      error('Exp.Pathways has %d entries per pathway, but the sequence has %d pulses.',size(sel_,2),N);
    end
  end
  for k = find(isRF)
    prev = 0;
    if k>1, prev = sel_(:,k-1); end
    if any(~isnan(sel_(:,k)) & ~isnan(prev) & sel_(:,k)~=prev)
      error('Exp.Pathways: RF pulses (pulse %d) cannot change the electron coherence order.',k);
    end
  end
  match = false(size(ctp,1),1);
  for r = 1:size(sel_,1)
    m_ = true(size(ctp,1),1);
    for k = 1:N
      if ~isnan(sel_(r,k))
        m_ = m_ & co(:,k)==sel_(r,k);
      end
    end
    match = match | m_;
  end
  ctp = ctp(match,:);
  co = co(match,:);
end
if isempty(ctp)
  error('No coherence transfer pathway is detectable (ending in -1) with Exp.Pathways.');
end

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
% Pathways that refocus at the detection point for all increments (echo
% pathways): at least one coherence order inversion, zero net coherence
% evolution for the initial delays and for each incremented dimension
function sel = echoPathways(co,t,Inc)
nP = size(co,1);
sel = false(nP,1);
for p = 1:nP
  c = co(p,:);
  cnz = c(c~=0);
  flipping = any(abs(diff(cnz))==2);
  t0 = c*t.';
  ok = flipping && abs(t0)<1e-12*max(1,max(t));
  for d = 1:max([Inc 0])
    ok = ok && sum(c(Inc==d))==0;
  end
  sel(p) = ok;
end
end

%-------------------------------------------------------------------------------
% Coherence transfer pathway parser: indices of the propagator blocks
% (1 = aa, 2 = bb, 3 = ab, 4 = ba) for the left and right side of each
% pulse and free evolution
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
