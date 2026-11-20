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
% nuclearRFFrame), with the rotating-wave approximation. For its error,
% rfDrive lists the largest field amplitude (mT) and the frequency band
% (MHz) of each RF pulse, one row per pulse and point.
rfEvents = find(cellfun(@(E) strcmp(E.type,'pulse') && isfield(E,'Channel') && strcmp(E.Channel,'rf'),Events));
useRFFrame = shiftedFrame && ~isempty(rfEvents);
rfDrive = zeros(0,3);
if useRFFrame
  logmsg(1,'  RF pulses in the nuclear RF frame (rotating-wave approximation)');
  for iEvent = rfEvents
    rfDrive(end+1,:) = [max(abs(Events{iEvent}.IQ(:))) Events{iEvent}.FrameBand]; %#ok<AGROW>
    iPulse = find(pulseEvents==iEvent);
    if isstruct(Vary) && isfield(Vary,'Pulses') && ~isempty(Vary.Pulses{iPulse})
      for k = 1:numel(Vary.Pulses{iPulse}.IQs)
        if ~isempty(Vary.Pulses{iPulse}.IQs{k})
          rfDrive(end+1,:) = [max(abs(Vary.Pulses{iPulse}.IQs{k}(:))) Vary.Pulses{iPulse}.FrameBands{k}]; %#ok<AGROW>
        end
      end
    end
  end
end

rawSignals = cell(1,nOrientations);
timeAxis = cell(1,nOrientations);
eventTimes = cell(1,nOrientations);
minWeight = ones(1,nOrientations);
minOverlap = ones(1,nOrientations);
rwaError = zeros(1,nOrientations);

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

  % Nuclear RF frame, the same for all RF pulses
  if useRFFrame
    if isempty(V)
      V = sf_eigenbasis(Ham);
    end
    [rfFrame,X] = nuclearRFFrame(Ham,V,-muxL,Sys.Spins,Sys.nElectrons);
    rwaError(iOrientation) = rfRWAError(rfFrame.E,rfFrame.f,X,rfDrive);
    if ~isempty(Relaxation_)
      rfFrame.Gamma = GammaEigen;
      rfFrame.equilibriumState = V'*Relaxation_.equilibriumState*V;
    end
    for iEvent = rfEvents
      Events_{iEvent}.rfFrame = rfFrame;
    end
  end

  [tAx, sig, ~, ~, ~, evT] = sf_propagate(Sigma_, Ham, Det, Events_, Relaxation_, Vary, Runs);
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
% Threshold on the relative error of the rotating-wave approximation for RF
% pulses (Bloch-Siegert shift relative to the nutation frequency, see
% rfRWAError). Measured deviation from lab-frame stepping (largest
% difference relative to the largest signal, single orientation) for the
% estimated error eps: Mims ENDOR of 1H at 350 mT with RF pi pulses of
% 2.5, 5, 10, 20 us: eps = 3.7e-3, 1.9e-3, 9.3e-4, 4.6e-4, deviation
% 9.6e-3, 3.1e-3, 1.3e-3, 9.4e-4; 14N with |A| > 2*nuL, 10 us: eps = 2.0e-3,
% deviation 1.6e-3; time-domain ENDOR with 1 us RF pi/2 pulses: eps =
% 4.2e-3, deviation 5.9e-3. The deviation is 0.8 to 2.6 times eps.
rwaTolerance = 0.01;
if useRFFrame
  logmsg(1,'  RF pulses: largest relative error of the rotating-wave approximation %0.2g',max(rwaError));
end
if max(rwaError)>rwaTolerance
  warning(['The rotating-wave approximation for RF pulses can be inaccurate ' ...
    '(estimated relative error %0.2g: RF frequency close to zero, or strong RF field). ' ...
    'Use Opt.SimFreq = 0 to simulate in the lab frame.'],max(rwaError));
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
% nuclearRFFrame  Eigenbasis and frame labels for RF pulses
%
%   [rfFrame,X] = nuclearRFFrame(H,V,X0,Spins,nElectrons)
%
%   RF pulses are propagated in a frame that rotates about the generator
%   G = V*diag(f)*V', which is diagonal in the eigenbasis V of H, so that H
%   commutes with G exactly and is not truncated (see rfFramePropagators in
%   sf_propagate). Each eigenstate k (column k of V, assigned to product
%   state k by sf_eigenbasis) gets the
%   label f(k) = sum over nuclei i of s_i*m_i(k), where m_i(k) is the m of
%   nucleus i in that product state, and s_i = +1 or -1 is the sense of
%   precession of nucleus i in the electron mS manifold of k: the sign of
%   tr(P_mS*Iz_i*H), the first-order z frequency of nucleus i in the
%   manifold (+1 if zero). The sense can differ between the manifolds, e.g.
%   if the hyperfine coupling is larger than twice the nuclear Zeeman
%   frequency. With these labels, the energy of the states rises with f
%   along the NMR transitions, so all of them rotate with the frame.
%
%   Input:
%     H           Hamiltonian (MHz) in the simulation frame, commuting with
%                 the total electron Sz
%     V           eigenvectors of H, from sf_eigenbasis
%     X0          RF excitation operator (MHz/mT), product basis: the field
%                 B1 (mT) along x(Lab) gives B1*X0
%     Spins       spin quantum numbers (Sys.Spins), electrons first
%     nElectrons  number of electron spins
%
%   Output:
%     rfFrame     structure with the fields
%       V           eigenvectors of H (columns), as input
%       E           eigenvalues (MHz), column
%       f           labels, column
%       Xp          elements of X with f(k)-f(l) = +1, zero elsewhere
%     X           RF excitation operator in the eigenbasis, V'*X0*V
function [rfFrame,X] = nuclearRFFrame(H,V,X0,Spins,nElectrons)
E = real(diag(V'*H*V));

% Electron manifold of each product state
mS = 0;
for iSpin = 1:nElectrons
  mS = mS + full(diag(sop(Spins,[iSpin,3])));
end
[~,~,manifold] = unique(mS);

% Labels: m of each nucleus, signed by its sense of precession in the
% manifold, tr(P_mS*Iz*H) = sum of m.*diag(H) over the manifold
Hdiag = real(diag(H));
f = zeros(size(H,1),1);
for iSpin = nElectrons+1:numel(Spins)
  m = full(diag(sop(Spins,[iSpin,3])));
  for iManifold = 1:max(manifold)
    inManifold = manifold==iManifold;
    s = sign(sum(m(inManifold).*Hdiag(inManifold)));
    if s==0, s = 1; end
    f(inManifold) = f(inManifold) + s*m(inManifold);
  end
end

X = V'*X0*V;
rfFrame.V = V;
rfFrame.E = E;
rfFrame.f = f;
rfFrame.Xp = X.*(abs(f-f.'-1)<1e-9);
end


%===============================================================================
% rfRWAError  Error of the rotating-wave approximation for RF pulses
%
%   err = rfRWAError(E,f,X,drive)
%
%   The RF field 2*Re(c*exp(2i*pi*F*t))*X0 is kept in the nuclear RF frame
%   only where it rotates with the frame (see nuclearRFFrame): the term
%   c*X(k,l) where f(k)-f(l) = -1, and its Hermitian conjugate. Relative to
%   H, a dropped term c*X(k,l) oscillates at E(k)-E(l)+F. Its leading
%   effect is a Bloch-Siegert shift, which, relative to the nutation
%   frequency, is the ratio of its amplitude |c*X(k,l)| to the frequency
%   distance |E(k)-E(l)+F|. The error is the largest ratio over all dropped
%   terms, with F anywhere in the frequency band of the pulse.
%
%   Input:
%     E      eigenvalues of H (MHz), column
%     f      labels of the eigenstates, column
%     X      RF excitation operator in the eigenbasis (MHz/mT)
%     drive  RF pulses, one row per pulse and point: largest field
%            amplitude (mT), lowest and highest frequency of the band (MHz)
%
%   Output:
%     err    largest relative error over all pulses (Inf if a dropped term
%            is resonant)
function err = rfRWAError(E,f,X,drive)
dropped = abs(f-f.'+1)>1e-9;
Xabs = abs(X).*dropped;
Xabs(Xabs<=1e-12*max(abs(X(:)))) = 0;
resonance = E.' - E;  % E(l)-E(k): frequency at which the term c*X(k,l) is static
err = 0;
for r = 1:size(drive,1)
  distance = max(max(drive(r,2)-resonance,resonance-drive(r,3)),0);
  ratio = drive(r,1)*Xabs./distance;
  ratio(Xabs==0) = 0;
  err = max(err,max(ratio(:)));
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
