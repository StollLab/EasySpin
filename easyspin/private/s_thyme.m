% thyme  Time domain evolution of density matrix
function [TimeArray, SignalArray, FinalStates, StateTrajectories, Events] = ...
  s_thyme(Sigma,Ham0,Det,Events,Relaxation,Vary)

if nargin==0, help(mfilename); return; end

if nargout>5, error('Too many output arguments.'); end
if nargin<4 || nargin>6, error('Wrong number of input arguments!'); end

nEvents = length(Events);

nDet = numel(Det.Ops);
Ham0 = Ham0*2*pi;
eigVectors = [];  % eigenbasis of Ham0, computed when needed

% Create some variables for bookkeeping
if ~isempty(Vary)
  nDimensions = numel(Vary.Points);
  DimensionIndices = ones(1,nDimensions);
  nPoints = prod(Vary.Points);
  IndirectDimensions = num2cell(Vary.Points); % This is used to store the output in structures that correspond to the dimensions
else
  nDimensions = 0;
  nPoints = 1;
  IndirectDimensions = {1}; % If no dimensions are provided, only one acquistion point exists
  DimensionIndices = 1;
end

n = size(Sigma,2);
FinalStates = zeros(IndirectDimensions{:},n,n);

%----------------------------------------------------------------------
% Detection: the signals are computed in the simulation frame for each
% detection operator D_k split into its parts D_k^(m) of electron
% coherence order m (one row of DetArray each, normalized by Tr(D_k'*D_k)),
% and combined in the detection frame by detectionFrameSignal
%----------------------------------------------------------------------
detOrders = unique(Det.M(:)).';
nOrders = numel(detOrders);
DetArray = zeros(nDet*nOrders,n^2);
for iDet = 1:nDet
  D = full(Det.Ops{iDet});
  normD = sum(abs(D(:)).^2);
  for iOrder = 1:nOrders
    Dm = D.*(Det.M==detOrders(iOrder));
    % row r with r*Sigma(:) = trace(Dm*Sigma)/normD
    DetArray((iDet-1)*nOrders+iOrder,:) = reshape(Dm.',1,n^2)/normD;
  end
end
nDetRows = size(DetArray,1);
StateTrajectories = [];
initialSigma = Sigma;
Resonator = false;

%----------------------------------------------------------------------
% Creates the pulse list, that is needed for reordering events
%----------------------------------------------------------------------
isPulse = zeros(1,nEvents);
nPulses = 0;
for iEvent = 1 : nEvents
  switch Events{iEvent}.type
    case 'pulse'
      isPulse(iEvent) = 1;
      nPulses = nPulses + 1;
      PulsePositions(nPulses) = iEvent; %#ok<AGROW>
      % If ringing is found, resonator is set to on
      if isfield(Events{iEvent},'Ringing')
        Resonator = true;
      end
  end
end

% Traces are stored in a numeric array, unless their lengths differ
% between acquisition points (see storage at the end of the loop)
StoreInCellArray = false;

% Initial event lengths, the indirect dimensions are applied to these
InitialEventLengths = zeros(1,nEvents);
for iEvent = 1 : nEvents
  InitialEventLengths(iEvent) = Events{iEvent}.t(end);
end

% Without indirect dimensions, the delay after a pulse is shortened by the
% (rounded) duration of the ringing of that pulse. s_sequencer has checked
% that the delays are long enough.
if nDimensions==0 && Resonator
  for iEvent = 1 : nEvents-1
    if strcmp(Events{iEvent}.type,'pulse')
      DecimalRound = 10;
      Ringing = round(Events{iEvent}.Ringing*(10^DecimalRound))/(10^DecimalRound);
      Events{iEvent+1}.t = Events{iEvent+1}.t - Ringing;
    end
  end
end

%----------------------------------------------------------------------
% Loop over number of Points = Product of the Points in each
% Dimension
%----------------------------------------------------------------------
for iPoints = 1 : nPoints
  
  Sigma = initialSigma;
  
  %--------------------------------------------------------------------
  % Set the pulse shapes and event lengths of the current acquisition
  % point
  %--------------------------------------------------------------------
  if ~isempty(Vary)
    % Apply the increments of the delays for the current acquisition
    EventLengths = InitialEventLengths;
    for iDimension = 1 : nDimensions
      if any(any(Vary.IncrementationTable{iDimension}))
        ModifiedEvents = find(Vary.IncrementationTable{iDimension}(:,DimensionIndices(iDimension)));
        for i = 1 : length(ModifiedEvents)
          EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable{iDimension}(ModifiedEvents(i),DimensionIndices(iDimension));
        end
      end
    end
    
    % Load the pulse shapes corresponding to the current acquisition
    % point (their propagators are updated below if the shape changed)
    for iPulse = 1 : nPulses
      if ~isempty(Vary.Pulses{iPulse})
        iEvent = PulsePositions(iPulse);
        Events{iEvent}.IQ = LoadfromArray(Vary.Pulses{iPulse}.IQs,DimensionIndices);
        Events{iEvent}.t = LoadfromArray(Vary.Pulses{iPulse}.ts,DimensionIndices);
        EventLengths(iEvent) = Events{iEvent}.t(end);
      end
    end
    
    % Reorder events and get new event lengths
    [Sequence, NewEventLengths] = s_reorder_events(EventLengths,isPulse);
    
    % Write new event durations of free evolution events to the event
    % structure
    for iEvent = 1 : nEvents
      NewPosition = find(Sequence == iEvent);
      if strcmp(Events{iEvent}.type,'free evolution') && NewEventLengths(NewPosition) ~= Events{iEvent}.t
        Events{iEvent}.t = NewEventLengths(NewPosition);
        Events{iEvent}.Propagation = [];
      end
    end
    
  else
    % If no Vary Structure is provided, events are processed linearly
    Sequence = 1 : nEvents;
  end
  
  % Start times of all events, from the start of the sequence
  tStart = eventStartTimes(Events,Sequence);
  
  %--------------------------------------------------------------------
  % Update the propagators of the pulses (free evolution propagators are
  % built during propagation)
  %--------------------------------------------------------------------
  % Each pulse gets the carrier phase at its start time (phase-continuous
  % carrier). Each step from t(k) to t(k+1) uses the field at its middle:
  % the two samples, moved by half a step at the center frequency of the
  % pulse (mean phase advance per step), are averaged.
  for iEvent = find(isPulse)
    IQ = Events{iEvent}.IQ.*exp(2i*pi*Events{iEvent}.FrameOffset*1e3*tStart(iEvent));
    if size(IQ,2)<2
      error('Pulse %d is shorter than the integration time step (Opt.IntTimeStep).',find(PulsePositions==iEvent));
    end
    halfStepPhaseRot = exp(1i*angle(sum(conj(IQ(:,1:end-1)).*IQ(:,2:end),'all'))/2);
    IQ = (IQ(:,1:end-1)*halfStepPhaseRot + IQ(:,2:end)*conj(halfStepPhaseRot))/2;
    Events{iEvent}.Propagation = pulsePropagators(Events{iEvent},IQ,Ham0,Relaxation);
  end
  
  %--------------------------------------------------------------------
  % Preparing for the actual propagation step
  %--------------------------------------------------------------------
  % Setting up some initial variables
  
  % bookkeping to keep track of when the first event is detected and
  % what should be stored, needed for proper arrangment of the output
  % later
  firstDetection = true;
  firstDensityMatrix = true;
  CurrentStateTrajectories = [];
  startTrace = 2;
  
  % Get the data points in all dimensions of the current acquistion,
  % will be used to store the output at the correct position in the
  % multidimensional arrays
  AcquisitionIndex = num2cell(DimensionIndices);
  
  %--------------------------------------------------------------------
  % Loop over events
  %--------------------------------------------------------------------
  for iEvent = Sequence
    
    currentEvent = Events{iEvent};
    
    %------------------------------------------------------------------
    % This computes the normalization values for the detection as well
    % as the first expectation value and the stores the state at the
    % beginning of this event in form of a density matrix if requested
    %------------------------------------------------------------------
    if currentEvent.Detection
      switch currentEvent.type
        case 'pulse'
          tvector = currentEvent.t;
        case 'free evolution'
          tvector = detectionTimes(currentEvent);
          if numel(tvector)>1
            currentEvent.TimeStep = tvector(2)-tvector(1);
          end
      end
      currentSignal = zeros(nDetRows,numel(tvector));
      currentSignal(:,1) = DetArray*Sigma(:);
      
    else
      switch currentEvent.type
        case 'pulse'
          tvector = currentEvent.t;
        case 'free evolution'
          tvector = [0 currentEvent.t];
          currentEvent.TimeStep = currentEvent.t;
      end
      
      currentSignal = [];
    end
    
    if currentEvent.StateTrajectories
      if currentEvent.Detection
        DensityMatrices  =cell(1,length(currentEvent.t));
      else
        DensityMatrices = cell(1,2);
      end
      DensityMatrices{1} = Sigma;
    else
      DensityMatrices = [];
    end
    %------------------------------------------------------------------
    
    switch currentEvent.type
      case 'pulse'
        %--------------------------------------------------------------
        % The time points of the pulse are the boundaries of the
        % propagation steps: the first point is the state before the
        % pulse, the last one the state at the end of the pulse.
        %--------------------------------------------------------------
        tvector = currentEvent.t;
        
        nPhaseSteps = size(currentEvent.PhaseCycle,1);
        
        %--------------------------------------------------------------
        % If phasecycling is requested for this event, a few variables
        % need to be defined
        %--------------------------------------------------------------
        if nPhaseSteps>1
          stateBeforePC = Sigma;
          SigmaPC = zeros(size(Sigma));
          PCnorm = sum(abs(currentEvent.PhaseCycle(:,2)));
          PCSignal = 0;
        end
        
        %--------------------------------------------------------------
        % Propagation Starts Here, and loops over all steps of the
        % phase cycle
        %--------------------------------------------------------------
        for iPhaseStep = 1 : nPhaseSteps
          
          if ~currentEvent.Relaxation
            nSteps = size(currentEvent.Propagation.Utotal,2);
          else
            nSteps = size(currentEvent.Propagation.Ltotal,2);
            SigmaVector = reshape(Sigma,n*n,1);
            if ~isempty(currentEvent.Propagation.PhaseTotal)
              phaseFactors = currentEvent.Propagation.PhaseTotal{iPhaseStep};
            end
          end
          
          %------------------------------------------------------------
          % Loop over the number of necessary propagation steps (1 if
          % detection during the pulse is not requested
          %------------------------------------------------------------
          for iWavePoint = 1 : nSteps
            if ~currentEvent.Relaxation
              U = currentEvent.Propagation.Utotal{iPhaseStep,iWavePoint};
              Sigma = U*Sigma*U';
            elseif currentEvent.Relaxation
              L = currentEvent.Propagation.Ltotal{iPhaseStep,iWavePoint};
              SigmaOffset = currentEvent.Propagation.SigmaOffsetTotal{iPhaseStep,iWavePoint};
              if isempty(currentEvent.Propagation.PhaseTotal)
                SigmaVector = L*SigmaVector + SigmaOffset;
              else
                % table Liouvillian rotated about the electron Sz (see pulsePropagators)
                % phase factors d(i)*conj(d(j)) of the vectorized state
                d = phaseFactors(:,iWavePoint);
                p = d.*d';
                p = p(:);
                SigmaVector = p.*(L*(conj(p).*SigmaVector) + SigmaOffset);
              end
              
              Sigma = reshape(SigmaVector,n,n);
            end
            
            %----------------------------------------------------------
            % Store expectation values and density matrix if requested
            %----------------------------------------------------------
            if currentEvent.Detection
              currentSignal(:,iWavePoint+1) = DetArray*Sigma(:);
            end
            
            if currentEvent.StateTrajectories
              DensityMatrices{iWavePoint+1} = Sigma;
            end
            %----------------------------------------------------------
          end
          
          %------------------------------------------------------------
          % Combine results from current phase cycle step with previous ones
          %------------------------------------------------------------
          if nPhaseSteps > 1
            PCweight = currentEvent.PhaseCycle(iPhaseStep,2)/PCnorm;
            SigmaPC = SigmaPC + PCweight*Sigma;
            if currentEvent.Detection
              PCSignal = PCSignal + PCweight*currentSignal;
            end
            if iPhaseStep ~= nPhaseSteps
              Sigma = stateBeforePC;
            end
          end
          %------------------------------------------------------------
          
        end
        
        %--------------------------------------------------------------
        % If phase cycling was requested, the results from combining
        % the individual signals and density matrices are used to
        % overwrite Sigma and currentSignal
        %--------------------------------------------------------------
        if nPhaseSteps > 1
          Sigma = SigmaPC;
          if currentEvent.Detection
            currentSignal = PCSignal;
          end
        end
        %--------------------------------------------------------------
        
      case 'free evolution'
        if currentEvent.Detection && ~currentEvent.Relaxation && ~currentEvent.StateTrajectories
          %--------------------------------------------------------------
          % Evolution in the eigenbasis of Ham0, for all detected time
          % points at once: element (k,l) of the state in the eigenbasis
          % evolves with exp(-1i*(E(k)-E(l))*t)
          %--------------------------------------------------------------
          if isempty(eigVectors)
            [eigVectors,E] = eig(full(Ham0+Ham0')/2,'vector');
            eigOmega = E-E.';
            eigOmega = eigOmega(:);
          end
          V = eigVectors;
          SigmaEig = V'*Sigma*V;
          % detection operators in the eigenbasis, using
          % sum(R.*Sigma) = sum((V.'*R*conj(V)).*SigmaEig)
          Coeff = zeros(nDetRows,n^2);
          for iRow = 1:nDetRows
            R = V.'*reshape(DetArray(iRow,:),n,n)*conj(V);
            Coeff(iRow,:) = R(:).'.*SigmaEig(:).';
          end
          chunkSize = max(1,floor(2^22/n^2));
          for k0 = 2:chunkSize:numel(tvector)
            k = k0:min(k0+chunkSize-1,numel(tvector));
            currentSignal(:,k) = Coeff*exp(-1i*eigOmega*tvector(k));
          end
          Sigma = V*(SigmaEig.*reshape(exp(-1i*eigOmega*tvector(end)),n,n))*V';

        else
          %--------------------------------------------------------------
          % Loads or computes the Propagator/Liouvillian if necessary
          %--------------------------------------------------------------
          if ~currentEvent.Relaxation
            if isfield(currentEvent.Propagation,'Utotal')
              U = currentEvent.Propagation.Utotal;
            else
              U = Propagator(Ham0,currentEvent.TimeStep);
              Events{iEvent}.Propagation.Utotal = U;
            end
          elseif currentEvent.Relaxation
            n = size(Sigma,1);
            SigmaVector = reshape(Sigma,n*n,1);
            if isfield(currentEvent.Propagation,'Ltotal')
              L = currentEvent.Propagation.Ltotal;
              SigmaOffset = currentEvent.Propagation.SigmaOffsetTotal;
            else
              Gamma = Relaxation.Gamma;
              equilibriumState = reshape(Relaxation.equilibriumState,n*n,1);
              [L, SigmaOffset] = Liouvillian(Ham0,Gamma,equilibriumState,currentEvent.TimeStep);
              Events{iEvent}.Propagation.Ltotal = L;
              Events{iEvent}.Propagation.SigmaOffsetTotal = SigmaOffset;
            end
          end
          %--------------------------------------------------------------

          %--------------------------------------------------------------
          % Loops over the time axis for propagation
          %--------------------------------------------------------------
          for itvector=2:length(tvector)
            if ~currentEvent.Relaxation
              Sigma = U*Sigma*U';
            elseif currentEvent.Relaxation
              SigmaVector = L*SigmaVector + SigmaOffset;
              Sigma = reshape(SigmaVector,n,n);
            end

            %------------------------------------------------------------
            % Stores expectation values and denisty matrices if requested
            %------------------------------------------------------------
            if currentEvent.Detection
              currentSignal(:,itvector) = DetArray*Sigma(:);
            end

            if currentEvent.StateTrajectories
              DensityMatrices{itvector} = Sigma;
            end
            %------------------------------------------------------------
          end
        end
    end
    
    %------------------------------------------------------------------
    % This combines the signals and time axes from all detected events
    % and prepares the output
    %------------------------------------------------------------------
    if ~isempty(currentSignal)
      currentSignal = detectionFrameSignal(currentSignal,tvector+tStart(iEvent),Det,detOrders);
      if firstDetection
        Signal = [];
        t = [];
        % store first point of the first signal
        Signal(:,1) = currentSignal(:,1);
        t(1,1) = tStart(iEvent);
        
        % now add all the others timepoints
        nSignal = size(currentSignal,2);
        endTrace = startTrace+nSignal-2;
        Signal(:,startTrace:endTrace) = currentSignal(:,2:end);
        t(1,startTrace:endTrace) = tvector(2:end)+tStart(iEvent);
        startTrace = endTrace+1;
        firstDetection = false;
        
      else
        % adding other signals. this slightly confusing index is
        % necessary in order to avoid double counting the last point of
        % the previously detected event and the first point of the
        % current signal
        nSignal = size(currentSignal,2);
        endTrace = startTrace+nSignal-2;
        Signal(:,startTrace:endTrace) = currentSignal(:,2:end);
        t(1,startTrace:endTrace) = tvector(2:end)+tStart(iEvent);
        startTrace = endTrace+1;
      end
    end
    
    % If StateTrajectories are required, they are stored in a temporary
    % cell array, which will later be stored to the output
    if currentEvent.StateTrajectories
      % Density matrices in the detection frame: an event stores either
      % all its time points or only its start and end
      if numel(DensityMatrices)==numel(tvector)
        tStates = tvector + tStart(iEvent);
      else
        tStates = [0 tvector(end)] + tStart(iEvent);
      end
      for iCell = 1 : length(DensityMatrices)
        DensityMatrices{iCell} = detectionFrameState(DensityMatrices{iCell},tStates(iCell),Det);
      end
      if firstDensityMatrix
        for iCell = 1 : length(DensityMatrices)
          CurrentStateTrajectories{iCell} = DensityMatrices{iCell}; %#ok<AGROW>
        end
        firstDensityMatrix = false;
        counter = length(DensityMatrices);
      else
        if currentEvent.Detection || iEvent == Sequence(end)
          iStart = counter;
          nElements = length(DensityMatrices);
          for iDensity = 2 : nElements
            CurrentStateTrajectories{iStart+iDensity-1} = DensityMatrices{iDensity}; %#ok<AGROW>
          end
          counter = counter + length(DensityMatrices) - 1;
        end
      end
    end
    %------------------------------------------------------------------
    
  end
  
  %--------------------------------------------------------------------
  % If time traces were stored they are now stored in a large numeric
  % array (or cell array if and of the detected events lengths changed)
  % that contains all the time traces according to the Vary structure
  %--------------------------------------------------------------------
  
  if ~firstDetection
    % If the current acquisition was the first, the arrays must first
    % be created
    if iPoints ==1
      if ~StoreInCellArray
        SignalSize = size(Signal);
        % Create an empty array for storage of all the signals
        SignalArray = zeros(IndirectDimensions{:},SignalSize(1),SignalSize(2));
        % Create an empty array for storage of all the time axes
        TimeArray = zeros(IndirectDimensions{:},SignalSize(2));
        
      else
        % Create cell arrays for output
        if isscalar(IndirectDimensions)
          % if no or only one Indirect Dimensions are requesteted, the
          % output structure is created here, to avoid creating square
          % cell arrays
          SignalArray = cell(1,IndirectDimensions{:});
          TimeArray = cell(1,IndirectDimensions{:});
        else
          SignalArray = cell(IndirectDimensions{:});
          TimeArray = cell(IndirectDimensions{:});
        end
      end
    end
    
    % Now the signals and their corresponding time axis can be stored.
    % The indexing uses AcquisitionIndex, which provides the current
    % position in the array
    if ~StoreInCellArray
      if length(t) > 1 && size(TimeArray,ndims(TimeArray)) ~= length(t)      
        % The trace of the current acquisition point has a different
        % length than the traces before, e.g. because a detected event is
        % incremented along an indirect dimension. From now on, traces are
        % stored in cell arrays, and the traces stored so far are moved
        % there. The number of points can also change if the total
        % detected time stays the same, since event lengths are rounded to
        % the time step, e.g. with dt = 0.2 µs:
        % 1st Acquisition: tau1 = 0.5 µs, tau2 = 0.5 µs, dt = 0.2 µs
        %                  total detection time = 0.5 + 0.5 = 1 µs
        %                  total number of data points = 2 + 2 = 4
        % 2nd Acquisition: tau1 = 0.4 µs, tau2 = 0.6 µs, dt = 0.2 µs
        %                  total detection time = 0.4 + 0.6 = 1 µs
        %                  total number of data points = 2 + 3 = 5
        StoreInCellArray = true;
        
        % Create cell arrays for output
        if isscalar(IndirectDimensions)
          % if no or only one Indirect Dimensions are requesteted, the
          % output structure is created here, to avoid creating square
          % cell arrays
          NewSignalArray = cell(1,IndirectDimensions{:});
          NewTimeArray = cell(1,IndirectDimensions{:});
        else
          NewSignalArray = cell(IndirectDimensions{:});
          NewTimeArray = cell(IndirectDimensions{:});
        end
        
        for jSignal = 1 : iPoints-1
          % get index for position of jSignal in the output cell array
          indexID = cell(length(IndirectDimensions),1);
          if isscalar(IndirectDimensions)
            indexID{1} = jSignal;
          else
            % acquisition points are visited with the last dimension
            % changing fastest
            [indexID{end:-1:1}] = ind2sub(fliplr(cell2mat(IndirectDimensions)),jSignal);
          end

          % assign jSignal to cell array, remove singleton dimensions
          NewSignalArray{indexID{:}} = squeeze(SignalArray(indexID{:},:,:));
          
          % If only one detection operator is used, the squeeze above also
          % rotates the signal vector and we have to reverse the rotation.
          if nDet==1
            NewSignalArray{indexID{:}} = NewSignalArray{indexID{:}}.';
          end
          
          % For the time axis, two cases need to be considered as well. If
          % there is only one indirect dimension, the time vector can be
          % taken as is
          if isscalar(IndirectDimensions)
            NewTimeArray{indexID{:}} = TimeArray(indexID{:},:);
          else
            % if there are more than one indirect dimensions, a squeeze has
            % to be applied and the vector needs to be rotated.
            NewTimeArray{indexID{:}} = squeeze(TimeArray(indexID{:},:)).';
          end
        end
        
        % store the signal and time axis of the most recent acquisition 
        % point (the one where the length of the transient changed)
        NewSignalArray{AcquisitionIndex{:}} = Signal;
        NewTimeArray{AcquisitionIndex{:}} = t;
        
        % Copy signal and time axes arrays over
        SignalArray = NewSignalArray;
        TimeArray = NewTimeArray;
                
      else 
        % Storing in vector array (transients have identical number of
        % points)
        TimeArray(AcquisitionIndex{:},:) = t;
        SignalArray(AcquisitionIndex{:},:,:) = Signal;
      end
    else
      % Storing in cell array (transients have different number of points)
      SignalArray{AcquisitionIndex{:}} = Signal;
      TimeArray{AcquisitionIndex{:}} = t;
    end
  end
  
  % Store the final state (in the detection frame) at its correct position
  tEnd = tStart(Sequence(end)) + Events{Sequence(end)}.t(end);
  FinalStates(AcquisitionIndex{:},:,:) = detectionFrameState(Sigma,tEnd,Det);
  
  % Store the cell array with state tractories (if any) of the current
  % acquisition point in its correct position in the output cellarray
  if ~isempty(CurrentStateTrajectories)
    StateTrajectories{AcquisitionIndex{:}} = CurrentStateTrajectories; %#ok<AGROW>
  end
  %--------------------------------------------------------------------
  % Incremeant the index for the Vary structure by 1
  %-------------------- ------------------------------------------------
  
  if ~isempty(Vary)
    for d = nDimensions:-1:1
      if DimensionIndices(d)<Vary.Points(d)
        DimensionIndices(d) = DimensionIndices(d)+1;
        break;
      else
        DimensionIndices(d) = 1;
      end
    end
  end
  
  %--------------------------------------------------------------------
end

%----------------------------------------------------------------------
% If detection was switched off during all events, the output for the
% time axis and signal is empty
%----------------------------------------------------------------------
if firstDetection
  SignalArray = [];
  TimeArray = [];
end
end


function tStart = eventStartTimes(Events,Sequence)
% Start times of all events in the order given by Sequence, with time
% counted from the start of the sequence. A pulse lasts until the last
% point of its time axis.
tStart = zeros(1,numel(Events));
t = 0;
for iEvent = Sequence
  tStart(iEvent) = t;
  t = t + Events{iEvent}.t(end);
end
end


function signal = detectionFrameSignal(rows,t,Det,detOrders)
% Signals in the detection frame, from the rows of DetArray*Sigma(:) (for
% each detection operator k, one row per electron coherence order m) at
% the times t (from the start of the sequence):
%   signal_k(t) = sum_m exp(-1i*m*(theta(t)+phi_k)) * row_{k,m}(t)
% with theta(t) = 2*pi*(DetFreq-SimFreq)*t and the detection phase phi_k
nOrders = numel(detOrders);
nDet = size(rows,1)/nOrders;
theta = 2*pi*Det.FrameOffset*1e3*t(:).';  % GHz*µs -> rad
signal = zeros(nDet,numel(t));
for iDet = 1:nDet
  iRows = (iDet-1)*nOrders + (1:nOrders);
  signal(iDet,:) = sum(exp(-1i*detOrders(:)*(theta+Det.Phase(iDet))).*rows(iRows,:),1);
end
end


function Sigma = detectionFrameState(Sigma,t,Det)
% Density matrix in the detection frame at time t (from the start of the
% sequence)
Sigma = Sigma.*exp(2i*pi*Det.FrameOffset*1e3*t*Det.M);
end


function t = detectionTimes(Event)
% Detected time points of a free evolution event, counted from its start:
% evenly spaced from 0 to Event.t, with steps of Event.DetTimeStep, or
% slightly shorter steps if the event length is not a multiple of it
if Event.t==0
  t = 0;
  return
end
nSteps = ceil(Event.t/Event.DetTimeStep-1e-9);
t = (0:nSteps)*(Event.t/nSteps);
t(end) = Event.t;
end


function Prop = pulsePropagators(Event,IQ,Ham0,Relaxation)
% Propagators (or Liouvillians, with relaxation) of a pulse event for each
% phase cycle step: one per time step, or, for a pulse without detection
% and relaxation, their product. IQ contains the field of each time step
% in the simulation frame (rows: phase cycle steps, columns: time steps),
% i.e. one value fewer per row than the time points in Event.t. The
% propagators are only rebuilt if IQ changed, and the lookup table of
% propagators is only rebuilt if its amplitudes changed.
Prop = Event.Propagation;
if isfield(Prop,'IQ') && isequal(Prop.IQ,IQ)
  return
end

nPhaseSteps = size(Event.PhaseCycle,1);

%------------------------------------------------------------
% convert the IQ wave form into a binary form, so that it
% is possible to use a propagator look up table
%------------------------------------------------------------
% The amplitude scale does not depend on the phase of the pulse, so that
% the lookup table can be reused for different start times
MaxWave = max(abs(IQ(:)));

if isfield(Event,'vertRes')
  % undocumented feature, you can provide your pulses with
  % their own vertical resolution e.g. for a rectangular
  % pulse it would suffice to have a vertRes of 3 (-1,0,1)
  if mod(Event.vertRes,2)
    vertRes = Event.vertRes+1;
  else
    vertRes = Event.vertRes;
  end
else
  vertRes = 1024;
end

% The pulse Hamiltonian is Re(IQ)*xOp + Im(IQ)*yOp (see s_propagationsetup),
% with a single operator if yOp is empty. There are three cases:
% - single operator: lookup table over Re(IQ), digitized in [-1 1]*MaxWave
% - rotated: in a shifted frame, the Hamiltonian is |IQ|*xOp rotated about
%   the total electron Sz by angle(IQ), and this rotation commutes with
%   Ham0. The propagator is obtained from a lookup table over |IQ|,
%   digitized in [0 1]*MaxWave, by multiplying its elements with
%   exp(-1i*angle(IQ)*M), M being the electron coherence order. With
%   relaxation, this requires that relaxation conserves the electron
%   coherence order and that the equilibrium state has no electron
%   coherences. Then the rotation is applied to the state during
%   propagation (PhaseTotal), and the table Liouvillians are not copied.
% - otherwise: propagator for each time step
twoOperators = ~isempty(Event.yOp);
rotated = twoOperators && ~isempty(Event.CoherenceOrder);

if Event.Relaxation % Set up for Liouville space simulations
  n = size(Ham0,1);
  Gamma = Relaxation.Gamma;
  equilibriumState = reshape(Relaxation.equilibriumState,n*n,1);
  if rotated
    Mvec = Event.CoherenceOrder(:);
    changesOrder = Mvec~=Mvec.';
    rotated = norm(Gamma(changesOrder)) <= 1e-10*norm(Gamma,'fro') && ...
      norm(equilibriumState(Mvec~=0)) <= 1e-10*norm(equilibriumState);
  end
end

if rotated
  amplitudes = 2*pi*MaxWave/vertRes*(0:vertRes);
  tableIndex = round(vertRes*abs(IQ)/MaxWave)+1;
  phase = angle(IQ);
elseif ~twoOperators
  % The factor 2 is due to the digitization of the wave
  amplitudes = 2*2*pi*MaxWave/vertRes*((0:vertRes)-vertRes/2);
  tableIndex = round(vertRes*(real(IQ)/MaxWave+1)/2)+1;
else
  realIQ = real(IQ);
  imagIQ = imag(IQ);
end
nWavePoints = size(IQ,2);
useTable = rotated || ~twoOperators;

%----------------------------------------------------------
% Build the Propagator table - increases speed by
% precomputing all propagators and later loading them,
% instead of calculating them at each way point (true as
% long as length(Wave) > vertical resolution)
%----------------------------------------------------------
if useTable && ~(isfield(Prop,'amplitudes') && isequal(Prop.amplitudes,amplitudes))
  if ~Event.Relaxation
    Prop.UTable = buildPropagators(Ham0,Event.xOp,Event.TimeStep,amplitudes);
  else
    [Prop.LTable, Prop.SigmaOffsetTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,Event.xOp,Event.TimeStep,amplitudes);
  end
  Prop.amplitudes = amplitudes;
end

% Set up structures for the actual propagators
if Event.Detection && ~Event.Relaxation
  Utotal = cell(nPhaseSteps,nWavePoints);
else
  if ~Event.Relaxation
    Utotal = cell(nPhaseSteps,1);
  else
    Ltotal = cell(nPhaseSteps,nWavePoints);
    SigmaOffsetTotal = Ltotal;
  end
end

% Local copies for the loop over the wave form (faster than field access)
relaxation = Event.Relaxation;
detection = Event.Detection;
if useTable
  if ~relaxation
    UTable = Prop.UTable;
  else
    LTable = Prop.LTable;
    SigmaOffsetTable = Prop.SigmaOffsetTable;
  end
end
if rotated
  M = Event.CoherenceOrder;
end
% Without detection and relaxation, only the product of the propagators
% is needed. The propagators of all time steps are taken from the table
% as a stack of matrices and multiplied with pagemtimes, in chunks to
% limit memory.
productOnly = useTable && ~relaxation && ~detection;
if productOnly
  n = size(Ham0,1);
  UTableStack = cat(3,UTable{:});
  chunkSize = max(2,floor(2^22/n^2));
end

for iPhaseStep = 1 : nPhaseSteps
  if rotated && ~relaxation
    % The rotation is diag(d)*U*diag(d)' with d = exp(-1i*angle(IQ)*mS),
    % where mS is taken relative to the first state (M(:,1)), since a
    % constant cancels
    D = exp(-1i*M(:,1)*phase(iPhaseStep,:));
  end
  if productOnly
    Uproduct = eye(n);
    for k0 = 1:chunkSize:nWavePoints
      k = k0:min(k0+chunkSize-1,nWavePoints);
      U = UTableStack(:,:,tableIndex(iPhaseStep,k));
      if rotated
        U = U.*reshape(D(:,k),n,1,[]).*reshape(conj(D(:,k)),1,n,[]);
      end
      Uproduct = productOfPropagators(U)*Uproduct;
    end
    Utotal{iPhaseStep,1} = Uproduct;
    continue
  end
  for iWavePoint = 1 : nWavePoints
    if useTable
      %----------------------------------------------------
      % Load from lookup table
      %----------------------------------------------------
      k = tableIndex(iPhaseStep,iWavePoint);
      if ~relaxation
        U = UTable{k};
        if rotated
          d = D(:,iWavePoint);
          U = (d.*U).*d';
        end
      else
        L = LTable{k};
        SigmaOffset = SigmaOffsetTable{k};
      end
    else
      %----------------------------------------------------
      % Compute propagator for this time step
      %----------------------------------------------------
      Ham1 = 2*pi*(realIQ(iPhaseStep,iWavePoint)*Event.xOp + imagIQ(iPhaseStep,iWavePoint)*Event.yOp);
      Ham =  Ham0+Ham1;
      if ~relaxation
        U = Propagator(Ham,Event.TimeStep);
      else
        [L, SigmaOffset] = Liouvillian(Ham,Gamma,equilibriumState,Event.TimeStep);
      end
    end

    if ~relaxation
      %------------------------------------------------------
      % Store Propagators for propation in one step
      % (additional speed bost if no detection) or stepwise
      %------------------------------------------------------
      if detection
        Utotal{iPhaseStep,iWavePoint} = U;
      else
        if iWavePoint == 1
          Uproduct = U;
        else
          Uproduct = U*Uproduct;
        end
      end
    else
      Ltotal{iPhaseStep,iWavePoint} = L;
      SigmaOffsetTotal{iPhaseStep,iWavePoint} = SigmaOffset;
    end
  end
  if ~relaxation && ~detection
    Utotal{iPhaseStep,1} = Uproduct;
  end
end

if ~Event.Relaxation
  Prop.Utotal = Utotal;
else
  Prop.Ltotal = Ltotal;
  Prop.SigmaOffsetTotal = SigmaOffsetTotal;
  if rotated
    % rotation d = exp(-1i*angle(IQ)*mS) for each time step, with mS
    % relative to the first state (M(:,1)), since a constant cancels
    Prop.PhaseTotal = cell(nPhaseSteps,1);
    for iPhaseStep = 1:nPhaseSteps
      Prop.PhaseTotal{iPhaseStep} = exp(-1i*M(:,1)*phase(iPhaseStep,:));
    end
  else
    Prop.PhaseTotal = [];
  end
end
Prop.IQ = IQ;
end


function P = productOfPropagators(U)
% Product U(:,:,end)*...*U(:,:,1) of a stack of propagators, obtained by
% repeatedly multiplying neighboring pairs
while size(U,3) > 1
  m = size(U,3);
  P = pagemtimes(U(:,:,2:2:m),U(:,:,1:2:m-1));
  if mod(m,2)
    P = cat(3,P,U(:,:,m));
  end
  U = P;
end
P = U;
end


function UTable = buildPropagators(Ham0,xOp,dt,amplitudes)

UTable = cell(1,numel(amplitudes));

for k = 1:numel(amplitudes)
  Ham = Ham0 + amplitudes(k)*xOp;
  UTable{k} = Propagator(Ham,dt);
end
end


function U = Propagator(Ham,dt)
U = expm(-1i*Ham*dt);
end


function [LTable, SigmaOffsetTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,xOp,dt,amplitudes)

SigmaOffsetTable = cell(1,numel(amplitudes));
LTable = cell(1,numel(amplitudes));

for k = 1:numel(amplitudes)
  Ham = Ham0 + amplitudes(k)*xOp;
  [LTable{k}, SigmaOffsetTable{k}] = Liouvillian(Ham,Gamma,equilibriumState,dt);
end
end


function [L, SigmaOffset] = Liouvillian(Ham,Gamma,equilibriumState,dt)

n = size(Ham,1);

HamSuOp = kron(eye(n,n),Ham)-kron(Ham.',eye(n,n));
L = -1i*HamSuOp-Gamma;

% Propagation over dt of d(sigma)/dt = L*sigma + Gamma*sigma_eq:
%   sigma(t+dt) = expm(L*dt)*sigma(t) + SigmaOffset
% with SigmaOffset = integral from 0 to dt of expm(L*s)*Gamma*sigma_eq ds.
% Both are obtained from the exponential of an augmented matrix, which
% works also if L is singular (e.g. due to trace conservation).
N = n^2;
A = zeros(N+1);
A(1:N,1:N) = L;
A(1:N,N+1) = Gamma*equilibriumState;
expA = expm(A*dt);
L = expA(1:N,1:N);
SigmaOffset = expA(1:N,N+1);
end


function LoadedElement = LoadfromArray(Array, ArrayIndex)
% This function loads an element from an Array (numeric or cell) for a
% given ArrayIndex (vector)
nPerDimension = size(Array);
IndexToLoad = ones(1,length(ArrayIndex));
if isscalar(ArrayIndex)
  if ArrayIndex <= nPerDimension(2)
    IndexToLoad = ArrayIndex;
  end
else
  for i = 1 : length(nPerDimension)
    if ArrayIndex(i) <= nPerDimension(i)
      IndexToLoad(i) = ArrayIndex(i);
    end
  end
end
IndexToLoad = num2cell(IndexToLoad);
if iscell(Array)
  LoadedElement = Array{IndexToLoad{:}};
else
  LoadedElement = Array(IndexToLoad{:});
end
end
