% sf_engine_fullspace  Full-space engine of saffron
%
%   [y,info] = sf_engine_fullspace(Sys,Exp,Opt,Orientations,nOrientations,averageOverChi)
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

function [y,info] = sf_engine_fullspace(Sys,Exp,Opt,Orientations,nOrientations,averageOverChi)

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

% Detection frame (validated in sf_propagationsetup)
if isfield(Exp,'DetFreq')
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

rawSignals = cell(1,nOrientations);
timeAxis = cell(1,nOrientations);
eventTimes = cell(1,nOrientations);
minWeight = ones(1,nOrientations);
minOverlap = ones(1,nOrientations);

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

  Relaxation_ = Relaxation;
  if ~isempty(Relaxation_)
    [Relaxation_.Gamma,minOverlap(iOrientation)] = sf_relaxationsuperoperator(Sys,Ham);
    if isempty(Relaxation_.equilibriumState)
      Relaxation_.equilibriumState = Sigma_;
    end
  end

  [tAx, sig, ~, ~, ~, evT] = sf_propagate(Sigma_, Ham, Det, Events, Relaxation_, Vary, Runs);
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

tensors = {'g', 'A', 'Q', 'D', 'ee', 'nn'};
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
