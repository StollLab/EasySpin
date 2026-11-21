% sf_engine_full  Full engine of saffron
%
%   [y,info] = sf_engine_full(Sys,Exp,Opt,Orientations,nOrientations,averageOverChi)
%
%   Propagates the density matrix of the full spin system in time (with
%   sf_propagate), for each orientation, and averages over orientations.
%
%   Input:
%     Sys             spin system (single isotopologue)
%     Exp             experiment, parsed by saffron; uses
%       Field           magnetic field (mT)
%       OriWeights      weights of the orientations
%       mwFreq, DetFreq, DetPhase  microwave and detection frame frequency
%                       (GHz), and detection phases
%       nPoints, DetWindow, DetIntegrate, Sequence
%     Opt             options, parsed by saffron, including
%       propEvents      events, with the wave forms of all settings of the
%                       pulses in the phase cycle, from sf_buildevents
%       propVary        indirect dimensions, from sf_buildevents
%       propCycleRuns   propagation runs of the phase cycle: settings of
%                       each pulse and signal weight, from phaseCycleRuns
%                       in saffron
%       propPathways    electron coherence orders of the selected pathways,
%                       from sf_pathwayorders in saffron
%       FrameShift      simulation frame frequency (GHz), from sf_buildevents
%       Relaxation      relaxation on/off, one value or one per event
%       relaxationFromSys  true if the user did not give Opt.Relaxation
%     Orientations    orientations (phi,theta,chi) in radians, one per row
%     nOrientations   number of orientations
%     averageOverChi  true if the orientations are averaged over the third
%                     Euler angle (powder); used for Sys.initState
%
%   Output:
%     y     signal, [indirect dimensions..., time, detection operator],
%           singleton dimensions removed
%     info  structure with fields
%       transientAxis  time axis of the transient (µs), empty for single
%                      points
%       EventTimes     [start end] of each element of Exp.Sequence (µs), one
%                      row per element, on the time axis of the transient

function [y,info] = sf_engine_full(Sys,Exp,Opt,Orientations,nOrientations,averageOverChi)

info = struct;
Events = Opt.propEvents;
Vary = Opt.propVary;
cycleRuns = Opt.propCycleRuns;
pathwayOrders = Opt.propPathways;
Opt = rmfield(Opt,{'propEvents','propVary','propCycleRuns','propPathways'});
nCycleRuns = numel(cycleRuns);

if ~isscalar(Opt.Relaxation) || Opt.relaxationFromSys
  for iEvent = 1:numel(Events)
    if isscalar(Opt.Relaxation)
      Events{iEvent}.Relaxation = Opt.Relaxation;
    else
      Events{iEvent}.Relaxation = Opt.Relaxation(iEvent);
    end
  end
end

logmsg(1,'-processing Sys structure------------------------------');

% Detection frame (validated in saffron)
if isfield(Exp,'DetFreq') && ~isempty(Exp.DetFreq)
  Opt.DetFreq = Exp.DetFreq;
elseif isfield(Exp,'mwFreq')
  Opt.DetFreq = Exp.mwFreq;
else
  Opt.DetFreq = 0;  % no microwave carrier: lab frame
end
if isfield(Exp,'DetPhase')
  Opt.DetPhase = Exp.DetPhase;
end

[Sys, Sigma, Det, Events, Relaxation] = sf_propagationsetup(Sys,Events,Opt);

% Simulation frame: shift by the frame frequency along the total electron Sz
shiftedFrame = Opt.FrameShift~=0;
if shiftedFrame
  logmsg(1,'  using secular Hamiltonian (shifted simulation frame)');
end
Sz_total = 0;
for iSpin = 1:Sys.nElectrons
  Sz_total = Sz_total + sop(Sys.Spins,[iSpin,3]);
end

% Coherence pathway selection: projection onto the electron coherence
% order after each pulse, one propagation run per selected pathway. Each
% propagation run is one run of the phase cycle with one pathway selection;
% sf_propagate computes the propagators once for all runs.
mS = full(diag(Sz_total));
Morder = mS - mS.';
pulseEvents = find(cellfun(@(E)strcmp(E.type,'pulse'),Events));
if isempty(pathwayOrders)
  pathwayOrders = nan(1,numel(pulseEvents));
end
if size(pathwayOrders,1)>1
  pathwayOrders = removeOverlaps(pathwayOrders,unique(round(Morder(:))));
end
nPathwayRuns = size(pathwayOrders,1);
nPulses = numel(pulseEvents);
Projections = cell(nPathwayRuns,nPulses);
for iPathway = 1:nPathwayRuns
  for k = 1:nPulses
    o = pathwayOrders(iPathway,k);
    if isinf(o)  % +1 or -1
      Projections{iPathway,k} = abs(abs(Morder)-1)<1e-10;
    elseif ~isnan(o)
      Projections{iPathway,k} = abs(Morder-o)<1e-10;
    end
  end
end
Runs = struct('rows',{},'factors',{},'Projection',{},'weight',{});
for iCycle = 1:nCycleRuns
  for iPathway = 1:nPathwayRuns
    Runs(end+1).rows = cycleRuns(iCycle).rows; %#ok<AGROW>
    Runs(end).factors = cycleRuns(iCycle).factors;
    Runs(end).Projection = Projections(iPathway,:);
    Runs(end).weight = cycleRuns(iCycle).weight;
  end
end

% User-provided initial state (Sys.initState), set up for each orientation
% (variables used in the parfor loop are defined in all cases)
userInitState = ~isempty(Sys.initState);
initStateBasis = '';
rho_M = [];
keep = true(Sys.nStates);
if userInitState
  if any(Opt.Relaxation)
    error('Sys.initState cannot be combined with relaxation in saffron.');
  end
  initStateBasis = Sys.initState{2};
  if strcmp(initStateBasis,'uncoupled')
    % molecular-frame state, rotated to the lab frame for each orientation
    rho_M = initstate(Sys);
  end
  % The orientational grid of a powder has only two angles (phi,theta), with
  % the third Euler angle chi set to zero. Rotating the molecule by chi
  % about the field is equivalent to shifting the phase of all microwave
  % pulses by -chi, and, with a secular Hamiltonian in a shifted frame, to
  % rotating the initial state and the detection operator about z by the
  % electron rotation exp(-1i*chi*Sz). An element rho(i,j) contributes with
  % exp(1i*p*chi), p = mS(i)-mS(j) the electron coherence order, so the chi
  % average keeps only the elements with p = 0. This is exact in a shifted
  % frame for microwave-only sequences.
  if averageOverChi
    keep = abs(mS-mS.')<1e-10;
  end
end

% RF pulses: their excitation operator is set up for each orientation. In
% a shifted frame, they are propagated in the nuclear RF frame (see
% sf_rfframe), with the rotating-wave approximation, if its error is below
% rwaTolerance, and are stepped otherwise (see rfSteppedPropagators in
% sf_propagate). For the error, rfDrive lists the largest field amplitude
% (mT) and the frequency band (MHz) of each RF pulse, one row per point,
% rfDriveEvent the RF pulse of each row, and rfDrivePoint the point (index
% into the wave forms of the pulse in Vary, 0 for the single wave form in
% Events).
% The choice is made for each point, so that e.g. an RF sweep is stepped
% only at the RF frequencies where the error is large.
% Threshold on the relative error of the rotating-wave approximation (the
% Bloch-Siegert shift relative to the nutation frequency). Measured
% deviation from lab-frame stepping (largest difference relative to the
% largest signal, single orientation) for the estimated error eps: Mims
% ENDOR of 1H at 350 mT with RF pi pulses of 2.5, 5, 10, 20 us: eps =
% 3.7e-3, 1.9e-3, 9.3e-4, 4.6e-4, deviation 9.6e-3, 3.1e-3, 1.3e-3, 9.4e-4;
% 14N with |A| > 2*nuL, 10 us: eps = 2.0e-3, deviation 1.6e-3; time-domain
% ENDOR with 1 us RF pi/2 pulses: eps = 4.2e-3, deviation 5.9e-3. The
% deviation is 0.8 to 2.6 times eps.
rwaTolerance = 0.01;
rfEvents = find(cellfun(@(E) strcmp(E.type,'pulse') && isfield(E,'Channel') && strcmp(E.Channel,'rf'),Events));
nRF = numel(rfEvents);
useRFFrame = shiftedFrame && nRF>0;
rfDrive = zeros(0,3);
rfDriveEvent = zeros(0,1);
rfDrivePoint = zeros(0,1);
rfPointsSize = cell(1,nRF);  % size of the array of wave forms in Vary, [] if none
if useRFFrame
  logmsg(1,'  RF pulses in the nuclear RF frame (rotating-wave approximation), stepped where its error exceeds %g',rwaTolerance);
  for k = 1:nRF
    iEvent = rfEvents(k);
    iPulse = find(pulseEvents==iEvent);
    if isstruct(Vary) && isfield(Vary,'Pulses') && ~isempty(Vary.Pulses{iPulse})
      % wave forms of the points (sf_propagate loads them for each point)
      rfPointsSize{k} = size(Vary.Pulses{iPulse}.IQs);
      for p = 1:numel(Vary.Pulses{iPulse}.IQs)
        if ~isempty(Vary.Pulses{iPulse}.IQs{p})
          rfDrive(end+1,:) = [max(abs(Vary.Pulses{iPulse}.IQs{p}(:))) Vary.Pulses{iPulse}.FrameBands{p}]; %#ok<AGROW>
          rfDrivePoint(end+1,1) = p; %#ok<AGROW>
        end
      end
    else
      rfDrive(end+1,:) = [max(abs(Events{iEvent}.IQ(:))) Events{iEvent}.FrameBand]; %#ok<AGROW>
      rfDrivePoint(end+1,1) = 0; %#ok<AGROW>
    end
    rfDriveEvent(end+1:size(rfDrive,1),1) = k;
  end
end

% Random RF phase (Exp.rfPhaseLocked = false): the RF carrier has a random
% phase in each shot, the same for all RF pulses of the shot. The signal
% is averaged over it with an internal cycle of nRFPhases equally spaced
% phases (see rfPhaseAverage). Within the rotating-wave approximation, a
% phase phi of the RF field rotates each RF propagator about G (see
% sf_rfframe), U -> R*U*R' with R = exp(-1i*phi*G). The state before an RF
% pulse and the detection operator propagated back to its end have parts
% with f orders between -K and K, K = max(f)-min(f) = 2*sum(I), which
% acquire the factors exp(1i*k*phi). Each RF pulse thus adds harmonics of
% phi up to order 2*K. RF pulses separated only by delays act as one,
% since the free evolution (and the projection onto electron coherence
% orders) commutes with G, so the rotations between them cancel. With
% nRFGroups such groups, the signal is a trigonometric polynomial of
% degree at most 2*nRFGroups*K, and its average over 2*nRFGroups*K+1
% equally spaced phases is exact. Without the rotating-wave approximation
% (stepped RF pulses, lab frame), the counter-rotating part of the field
% adds higher harmonics, and the number of phases is doubled until the
% result converges.
rfRandomPhase = nRF>0 && ~(isfield(Exp,'rfPhaseLocked') && Exp.rfPhaseLocked);
if rfRandomPhase
  % groups of RF pulses without microwave pulses between them, in the order
  % of the sequence; if a delay can become negative along an indirect
  % dimension, pulses can change their order, and each RF pulse counts
  isRFPulse = ismember(pulseEvents,rfEvents);
  nRFGroups = sum(diff([false isRFPulse])==1);
  if isstruct(Vary)
    minLength = cellfun(@(E) E.t(end),Events);
    for d = 1:numel(Vary.IncrementationTable)
      increments = Vary.IncrementationTable{d};  % events x points, empty if none
      if ~isempty(increments)
        minLength = minLength + min(0,min(increments,[],2)).';
      end
    end
    if any(minLength<0), nRFGroups = nRF; end
  end
  nRFPhases = 2*nRFGroups*sum(2*Sys.Spins(Sys.nElectrons+1:end)) + 1;
  logmsg(1,'  random RF phase: average over %d phases',nRFPhases);
else
  nRFPhases = 1;
end

rawSignals = cell(1,nOrientations);
timeAxis = cell(1,nOrientations);
eventTimes = cell(1,nOrientations);
minWeight = ones(1,nOrientations);
minOverlap = ones(1,nOrientations);
rwaError = zeros(1,nOrientations);
nSteppedRF = zeros(1,nOrientations);
nPhasesUsed = ones(1,nOrientations);
phaseConverged = true(1,nOrientations);

logmsg(1,'-starting orientation loop-----------------------------');
logmsg(1,'  propagating %d orientations, %d phase cycle run(s), %d pathway run(s)',...
  nOrientations,nCycleRuns,nPathwayRuns);
Field = Exp.Field;

% A parallel pool only pays off for several orientations (starting one
% takes seconds), so a single orientation runs without one.
if nOrientations>1, maxWorkers = Inf; else, maxWorkers = 0; end

parfor (iOrientation = 1 : nOrientations, maxWorkers)

  Sys_L = rotatesystem(Sys,Orientations(iOrientation,:)); % mol->lab frame

  Ham = ham(Sys_L,Field*[0 0 1]);  % lab frame

  % RF excitation operator (MHz/mT): nuclear Zeeman interaction (with
  % gnscale and chemical shielding) with an RF field along x(Lab), H1 =
  % 2*Re(IQ)*X0 with the lab-frame RF waveform IQ (mT) and X0 = -mux. It
  % has electron coherence order 0, so it is the same in shifted frames.
  % The hyperfine enhancement via the electron Zeeman interaction with the
  % RF field is not included.
  Events_ = Events;
  muxL = [];
  if ~isempty(rfEvents)
    [muxL,~,~] = ham_nz(Sys_L);
    for iEvent = rfEvents
      Events_{iEvent}.xOp = -2*muxL;
    end
  end
  if shiftedFrame
    [Ham,mSweights] = sf_secularham(Ham,Sys.Spins,Sys.nElectrons);  % high-field approximation
    minWeight(iOrientation) = min(mSweights);
  end

  % Initial state in the lab frame
  Sigma_ = Sigma;
  if userInitState
    switch initStateBasis
      case 'uncoupled'
        U = spinrotation(Sys.Spins,Orientations(iOrientation,:));
        Sigma_ = U'*rho_M*U;
      case 'eigen'
        Sigma_ = initstate(Sys_L,Ham);
      case 'zerofield'
        Sigma_ = initstate(Sys_L);
    end
    Sigma_(~keep) = 0;
  end

  Ham = Ham - Opt.FrameShift*1e3*Sz_total;  % simulation frame (MHz)

  % Relaxation superoperator, in the product basis and in the eigenbasis V
  % of Ham (the latter for the nuclear RF frame)
  Relaxation_ = Relaxation;
  V = [];
  GammaEigen = [];
  if ~isempty(Relaxation_)
    [Relaxation_.Gamma,minOverlap(iOrientation),GammaEigen,V] = sf_relaxationsuperoperator(Sys,Ham);
    if isempty(Relaxation_.equilibriumState)
      Relaxation_.equilibriumState = Sigma_;
    end
  end

  % Nuclear RF frame, the same for all RF pulses. An RF pulse is stepped in
  % the shifted microwave frame instead at the points where its
  % rotating-wave error exceeds the tolerance. This is exact, since the RF
  % operator has electron coherence order 0 and is the same in that frame.
  stepped = false(size(rfDrive,1),1);
  if useRFFrame
    if isempty(V)
      V = sf_eigenbasis(Ham);
    end
    [rfFrame,err] = sf_rfframe(Ham,V,-muxL,Sys.Spins,Sys.nElectrons,rfDrive);
    stepped = err>rwaTolerance;
    rwaError(iOrientation) = max(err);
    nSteppedRF(iOrientation) = sum(stepped);
    if ~isempty(Relaxation_)
      rfFrame.Gamma = GammaEigen;
      rfFrame.equilibriumState = V'*Relaxation_.equilibriumState*V;
    end
    for k = 1:nRF
      rows = rfDriveEvent==k;
      if isempty(rfPointsSize{k})
        steppedPoints = stepped(rows);  % single wave form
      else
        % one value per wave form in Vary, read by sf_propagate per point
        steppedPoints = false(rfPointsSize{k});
        steppedPoints(rfDrivePoint(rows)) = stepped(rows);
      end
      Events_{rfEvents(k)}.rfFrame = rfFrame;
      Events_{rfEvents(k)}.rfStepped = steppedPoints;
    end
  end

  % Average over the random RF phase: exact with nRFPhases phases if all RF
  % pulses use the rotating-wave approximation, otherwise until converged
  adaptive = ~useRFFrame || any(stepped);
  [tAx,sig,evT,nPhasesUsed(iOrientation),phaseConverged(iOrientation)] = ...
    rfPhaseAverage(Sigma_,Ham,Det,Events_,Relaxation_,Vary,Runs,rfEvents,pulseEvents,nRFPhases,adaptive);
  timeAxis{iOrientation} = tAx;
  rawSignals{iOrientation} = sig;
  eventTimes{iOrientation} = evT;
end

if min(minWeight)<0.9
  warning(['The electron spin states are strongly mixed (smallest weight in an ' ...
    'mS manifold: %0.2f). The high-field approximation used in the shifted ' ...
    'simulation frame can be inaccurate. Use Opt.SimFreq = 0 to simulate in ' ...
    'the lab frame.'],min(minWeight));
end
if useRFFrame
  logmsg(1,'  RF pulses: largest relative error of the rotating-wave approximation %0.2g',max(rwaError));
  if any(nSteppedRF)
    logmsg(1,'  RF pulses stepped (error above %g) at some points in %d of %d orientations',...
      rwaTolerance,sum(nSteppedRF>0),nOrientations);
  end
end
if rfRandomPhase && any(nPhasesUsed>nRFPhases)
  logmsg(1,'  random RF phase: up to %d phases (counter-rotating RF field)',max(nPhasesUsed));
end
if ~all(phaseConverged)
  warning(['The average over the random RF phase did not converge with %d phases ' ...
    '(very strong RF field at low frequencies). The result can be inaccurate.'],max(nPhasesUsed));
end
relaxationMatrix = ~isscalar(Sys.T1) || ~isscalar(Sys.T2);
if relaxationMatrix && min(minOverlap)<0.5
  warning(['The eigenstates are strongly mixed (smallest overlap with a product ' ...
    'state: %0.2f). The assignment of the relaxation times in Sys.T1 and Sys.T2 ' ...
    'to transitions is ambiguous.'],min(minOverlap));
end

timeAxis = timeAxis{1};

% Integrate signal over orientations
Signal = 0;
for iOrientation = 1 : nOrientations
  Signal = Signal + rawSignals{iOrientation}*Exp.OriWeights(iOrientation);
end

% Arrange as [indirect dimensions..., time, detection operator]
if isfield(Exp,'nPoints') && ~isempty(Exp.nPoints)
  nPts = Exp.nPoints(:).';
else
  nPts = 1;
end
nDet = numel(Det.Ops);
nTime = numel(Signal)/prod(nPts)/nDet;
Signal = reshape(Signal,[prod(nPts) nDet nTime]);
Signal = permute(Signal,[1 3 2]);  % [points, time, det]

% detected time points of the first acquisition (from the start of the sequence)
if iscell(timeAxis)
  tDetected = timeAxis{1};
else
  tDetected = reshape(timeAxis,[],size(timeAxis,ndims(timeAxis)));
  tDetected = tDetected(1,:);
end

integrateEchoes = isfield(Exp,'DetIntegrate') && Exp.DetIntegrate && nTime>1;
if integrateEchoes
  % integral over the detection window (trapezoidal rule)
  Signal = trapz(tDetected,Signal,2);
  nTime = 1;
end

y = reshape(Signal,[nPts nTime nDet 1]);
y = squeezeSignal(y);

% Transient axis
if nTime>1
  if isnumeric(Exp.DetWindow)
    info.transientAxis = Exp.DetWindow(1) + (tDetected - tDetected(1));
  else
    info.transientAxis = tDetected;  % whole sequence, from its start
  end
else
  info.transientAxis = [];
end

% Start and end of each element of Exp.Sequence, on the time axis of the
% transient: from the start of the sequence for DetWindow = 'sequence',
% from the end of the sequence otherwise. Events appended for the
% detection window are left out; the last of them starts at DetWindow(1).
nElements = numel(Exp.Sequence);
evT = eventTimes{1};
if isnumeric(Exp.DetWindow)
  tEndSequence = evT(end,1) - Exp.DetWindow(1);
  evT = evT(1:nElements,:) - tEndSequence;
  if ~isstruct(Exp.Sequence{end})
    evT(end,2) = 0;  % the last delay, shortened or extended by DetWindow(1)
  end
else
  evT = evT(1:nElements,:);
end
info.EventTimes = evT;

end


%===============================================================================
% removeOverlaps  Replace overlapping pathway selections by separate pathways
%
%   orders = removeOverlaps(orders,allOrders)
%
%   Each row of orders is propagated separately and the signals are added, so
%   a pathway selected by two rows would count twice. If any rows overlap,
%   the selections are expanded into fully specified pathways, each listed
%   once. Selections without overlap are returned unchanged.
%
%   Input:
%     orders     electron coherence order after each pulse (NaN = any, Inf =
%                +1 or -1), one row per selection and one column per pulse
%     allOrders  all electron coherence orders of the spin system, column
%
%   Output:
%     orders     selections without overlap
function orders = removeOverlaps(orders,allOrders)
expanded = zeros(0,size(orders,2));
for r = 1:size(orders,1)
  rowPathways = zeros(1,0);
  for k = 1:size(orders,2)
    o = orders(r,k);
    if isnan(o)
      values = allOrders;
    elseif isinf(o)
      values = [-1; +1];
    else
      values = o;
    end
    rowPathways = [repelem(rowPathways,numel(values),1) repmat(values,size(rowPathways,1),1)];
  end
  expanded = [expanded; rowPathways]; %#ok<AGROW>
end
distinct = unique(expanded,'rows');
if size(distinct,1)<size(expanded,1)
  orders = distinct;
end
end


%===============================================================================
% rfPhaseAverage  Signal averaged over the random phase of the RF carrier
%
%   [t,sig,evT,nPhases,converged] = rfPhaseAverage(Sigma,Ham,Det,Events,Relaxation,Vary,Runs,rfEvents,pulseEvents,nPhases,adaptive)
%
%   Propagates with the RF carrier phases 2*pi*j/nPhases, j = 0..nPhases-1,
%   added to all RF pulses (see rfPhaseRuns), and averages the signals.
%   With adaptive = true, the number of phases is doubled, adding the
%   phases halfway between the earlier ones, until the averages over
%   nPhases and 2*nPhases phases differ by less than phaseTolerance
%   relative to the largest signal. The average over the doubled number of
%   phases is returned.
%
%   Input:
%     Sigma, Ham, Det, Events, Relaxation, Vary, Runs  inputs of sf_propagate
%     rfEvents     indices of the RF pulse events in Events
%     pulseEvents  indices of the pulse events in Events
%     nPhases      number of phases (1: single phase, no average)
%     adaptive     true to double nPhases until the average converges
%
%   Output:
%     t, sig, evT  time axis, signal and event times, from sf_propagate
%     nPhases      number of phases used
%     converged    false if the adaptive average did not converge
function [t,sig,evT,nPhases,converged] = rfPhaseAverage(Sigma,Ham,Det,Events,Relaxation,Vary,Runs,rfEvents,pulseEvents,nPhases,adaptive)
converged = true;
if nPhases==1
  [t,sig,~,~,~,evT] = sf_propagate(Sigma,Ham,Det,Events,Relaxation,Vary,Runs);
  return
end
phases = 2*pi*(0:nPhases-1)/nPhases;
[E,V,R] = rfPhaseRuns(Events,Vary,Runs,rfEvents,pulseEvents,phases);
[t,sig,~,~,Eout,evT] = sf_propagate(Sigma,Ham,Det,E,Relaxation,V,R);
if ~adaptive || isempty(sig)
  return
end
phaseTolerance = 1e-4;
maxPhases = 64*nPhases;
while true
  % Phases halfway between the earlier ones. The pulse propagators are
  % kept (they are rebuilt only where the wave form changed); other fields
  % of the events are reset, since sf_propagate changes event lengths.
  for iEvent = pulseEvents
    Events{iEvent}.Propagation = Eout{iEvent}.Propagation;
  end
  [E,V,R] = rfPhaseRuns(Events,Vary,Runs,rfEvents,pulseEvents,phases+pi/nPhases);
  [~,sigNew,~,~,Eout] = sf_propagate(Sigma,Ham,Det,E,Relaxation,V,R);
  sigDoubled = (sig+sigNew)/2;
  change = max(abs(sigDoubled(:)-sig(:)));
  sig = sigDoubled;
  phases = 2*pi*(0:2*nPhases-1)/(2*nPhases);
  nPhases = 2*nPhases;
  if change<=phaseTolerance*max(abs(sig(:)))
    break
  end
  if nPhases>=maxPhases
    converged = false;
    break
  end
end
end


%===============================================================================
% rfPhaseRuns  Events and propagation runs for several RF carrier phases
%
%   [Events,Vary,Runs] = rfPhaseRuns(Events,Vary,Runs,rfEvents,pulseEvents,phases)
%
%   The settings of each RF pulse (rows of IQ and PhaseCycle, also for the
%   points of the indirect dimensions in Vary) are repeated for each phase
%   phi in phases, with the wave forms multiplied by exp(1i*phi). The
%   signal is averaged over the phases, with weights 1/numel(phases):
%   - one RF pulse: the states after the pulse are averaged over its
%     settings for all phases within each run (see sf_propagate)
%   - several RF pulses: one run per phase, in which all RF pulses use the
%     settings of that phase, since they share the phase of the carrier
%
%   Input:
%     Events       events, from sf_buildevents
%     Vary         indirect dimensions, from sf_buildevents, or []
%     Runs         propagation runs (see sf_propagate)
%     rfEvents     indices of the RF pulse events in Events
%     pulseEvents  indices of the pulse events in Events
%     phases       RF carrier phases (radians)
%
%   Output:
%     Events, Vary, Runs  with the settings and runs of all phases
function [Events,Vary,Runs] = rfPhaseRuns(Events,Vary,Runs,rfEvents,pulseEvents,phases)
nPhases = numel(phases);
phaseFactors = exp(1i*phases(:));
nSettings = zeros(1,numel(rfEvents));
iPulses = zeros(1,numel(rfEvents));
for k = 1:numel(rfEvents)
  iEvent = rfEvents(k);
  iPulse = find(pulseEvents==iEvent);
  iPulses(k) = iPulse;
  PC = Events{iEvent}.PhaseCycle;
  nSettings(k) = size(PC,1);
  Events{iEvent}.PhaseCycle = repmat(PC,nPhases,1) + repelem(phases(:),nSettings(k),size(PC,2));
  Events{iEvent}.IQ = kron(phaseFactors,Events{iEvent}.IQ);
  if isstruct(Vary) && isfield(Vary,'Pulses') && ~isempty(Vary.Pulses{iPulse})
    Vary.Pulses{iPulse}.IQs = cellfun(@(IQ) kron(phaseFactors,IQ),Vary.Pulses{iPulse}.IQs,'UniformOutput',false);
  end
end

if isscalar(rfEvents)
  p = iPulses;
  for r = 1:numel(Runs)
    Runs(r).rows{p} = reshape(Runs(r).rows{p}(:) + (0:nPhases-1)*nSettings,[],1);
    Runs(r).factors{p} = repmat(Runs(r).factors{p}(:),nPhases,1)/nPhases;
  end
else
  RunsPhases = repmat(Runs(:),1,nPhases);
  for r = 1:numel(Runs)
    for j = 1:nPhases
      for k = 1:numel(rfEvents)
        p = iPulses(k);
        RunsPhases(r,j).rows{p} = Runs(r).rows{p} + (j-1)*nSettings(k);
      end
      RunsPhases(r,j).weight = Runs(r).weight/nPhases;
    end
  end
  Runs = reshape(RunsPhases.',1,[]);
end
end


%===============================================================================
% squeezeSignal  Remove singleton dimensions of the signal
%
%   y = squeezeSignal(y)
%
%   Input:
%     y  signal array
%
%   Output:
%     y  signal array without singleton dimensions; 1D results are row
%        vectors
function y = squeezeSignal(y)
sz = size(y);
sz(sz==1) = [];
if isempty(sz)
  y = y(1);
elseif isscalar(sz)
  y = reshape(y,1,sz);
else
  y = reshape(y,sz);
end
end

%===============================================================================
% rotatesystem  Transform a spin system from the molecular to the lab frame
%
%   Sys_L = rotatesystem(Sys_M,angles_M2L)
%
%   Input:
%     Sys_M       spin system in the molecular frame
%     angles_M2L  Euler angles [phi theta chi] (radians) of the
%                 transformation from the molecular frame to the lab frame
%
%   Output:
%     Sys_L       spin system with all tensor frames (gFrame, AFrame, ...)
%                 relative to the lab frame
function Sys_L = rotatesystem(Sys_M,angles_M2L)

Sys_L = Sys_M;

R_M2L = erot(angles_M2L);  % molecular frame -> lab frame
R_L2M = R_M2L.'; % lab frame -> molecular frame

tensors = {'g', 'A', 'Q', 'D', 'ee', 'nn', 'sigma'};
tensors2rotate = tensors(isfield(Sys_M,tensors));

for t = 1:numel(tensors2rotate)
  tFrame_ = [tensors2rotate{t} 'Frame'];
  ang_M2T = Sys_M.(tFrame_);
  ang_L2T = zeros(size(ang_M2T));
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
% spinrotation  Spin rotation operator for the mol->lab frame transformation
%
%   U = spinrotation(Spins,angles_M2L)
%
%   Input:
%     Spins       spin quantum numbers of all spins (Sys.Spins)
%     angles_M2L  Euler angles [phi theta chi] (radians) of the
%                 transformation from the molecular frame to the lab frame
%
%   Output:
%     U           rotation operator such that H_L = U'*H_M*U and
%                 rho_L = U'*rho_M*U, consistent with rotatesystem
function U = spinrotation(Spins,angles_M2L)
U = 1;
for S = Spins
  U = kron(U,wignerd(S,angles_M2L(1),angles_M2L(2),angles_M2L(3)));
end
end
