% thyme  Time domain evolution of density matrix
function [TimeArray, SignalArray, FinalStates, StateTrajectories, Events] = s_thyme(Sigma,Ham0,Det,Events,Relaxation,Vary)

if (nargin==0), help(mfilename); return; end

if (nargout>5), error('Too many output arguments.'); end
if (nargin<4) || (nargin>6), error('Wrong number of input arguments!'); end


nEvents = length(Events);


nDet = numel(Det);
normsDet = zeros(1,nDet);
Ham0 = Ham0*2*pi;

% Create some variables for bookkeeping
if ~isempty(Vary)
  nDimensions = numel(Vary.Points);
  DimensionIndices = ones(1,nDimensions);
  nPoints = prod(Vary.Points);
  IndirectDimensions = num2cell(Vary.Points); % This is used to store the out put in structures that correspond to the dimensions
else
  nDimensions = 0;
  nPoints = 1;
  IndirectDimensions = {1}; % If no dimensions are provided, only one acquistion point exists
  DimensionIndices = 1;
end

n = size(Sigma,2);
FinalStates = zeros(IndirectDimensions{:},n,n);
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

%----------------------------------------------------------------------
% This checks if traces for all datapoints have the save length, if
% not, a cell array is used for storage of traces
%----------------------------------------------------------------------
StoreInCellArray = false;

% Gets initial event lengths
InitialEventLengths = zeros(1,nEvents);
for iEvent = 1 : nEvents
  InitialEventLengths(iEvent) = Events{iEvent}.t(end);
end

if nDimensions > 0
  
  
  % Loop over all data points
  for iPoint = 1 : nPoints
    
    EventLengths = InitialEventLengths;
    
    % DetectionTime is used to measure the total length of detected
    % events
    DetectionTime = 0;
    
    % First all delays are set according to the incrementation tables
    % for each dimension
    for iDimension = 1 : nDimensions
      if any(any(Vary.IncrementationTable{iDimension}))
        ModifiedEvents = find(Vary.IncrementationTable{iDimension}(:,DimensionIndices(iDimension)));
        for i = 1 : length(ModifiedEvents)
          EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable{iDimension}(ModifiedEvents(i),DimensionIndices(iDimension));
        end
      end
    end
    
    % Loops over the pulses, and loads the pulse length for the current
    % acquisition point
    for iPulse = 1 : nPulses
      if ~isempty(Vary.Pulses{iPulse})
        iEvent = PulsePositions(iPulse);
        t = LoadfromArray(Vary.Pulses{iPulse}.ts,DimensionIndices);
        if Resonator
          % if resonator is present the original length of the pulse is
          % calculated (and rounded)
          DecimalRound = 10;
          Ringing = LoadfromArray(Vary.Pulses{iPulse}.Ringing,DimensionIndices);
          Ringing = round(Ringing*(10^DecimalRound))/(10^DecimalRound);
          EventLengths(iEvent) = t(end)-Ringing;
        else
          % if not, the length can be obtained from the last element of
          % the time axis
          EventLengths(iEvent) = t(end);
        end
      end
    end
    
    % Reorder Events and adjust event lengths
    [Sequence, NewEventLengths] = s_reorder_events(EventLengths,isPulse);
    
    if Resonator
      % Loop over the pulses
      for iPulse = 1 : nPulses
        % Get the initial event index of the current pulse
        ThisEventIndex = PulsePositions(iPulse);
        % Use the inital event index to find it in the position after
        % reordering
        ThisEvent = find(Sequence == PulsePositions(iPulse));
        if ~isempty(Vary.Pulses{iPulse})
          % if any parameter of the pulse is varied, the length of the
          % ringing is taken from the vary structure
          Ringing = LoadfromArray(Vary.Pulses{iPulse}.Ringing,DimensionIndices);
        else
          % if the pulse is not in the vary structure, the length of
          % the ringing is taken from the Event structure
          Ringing = Events{Sequence(ThisEventIndex)}.Ringing;
        end
        % Rounding of the Ringing length
        DecimalRound = 10;
        Ringing = round(Ringing*(10^DecimalRound))/(10^DecimalRound);
        % Length of the pulse is being updated, to make sure that the
        % counter for detection periods is correct
        NewEventLengths(ThisEvent) = NewEventLengths(ThisEvent) + Ringing;
        if ThisEvent < nEvents
          % if the current pulse is not the last pulse in the new
          % sequence, the succeeding event is checked to make sure it
          % is not a pulse
          FollowingEvent = Sequence(ThisEvent+1);
          if strcmp(Events{FollowingEvent}.type,'pulse')
            error('When using a resonator, pulses need to be separated by inter pulse delays to accomodate for ringing from the resonator.')
          end
          % And then the free evolution event is shortened by the
          % length of Ringing
          NewEventLengths(FollowingEvent) = NewEventLengths(FollowingEvent) - Ringing;
          % And the new length is checked again
          if NewEventLengths(FollowingEvent) < 0
            Msg = ['The delay ' num2str(FollowingEvent) ' is too short to accomodate for ringing of the preceding pulse.'];
            error(Msg);
          end
        end
      end
    end
    
    % Count the time during detected events
    for iSequence = 1 : nEvents
      if Events{Sequence(iSequence)}.Detection
        DetectionTime = DetectionTime + NewEventLengths(iSequence);
      end
    end
    
    % Get a reference detection time for the first data point
    if iPoint == 1
      ReferenceDetection = DetectionTime;
    else
      % Or compare the sum of the detected time intervals to the
      % detection time of the first data point. If they don't match,
      % signals are stored in a cell array
      if abs(ReferenceDetection-DetectionTime) > 1e-12
        StoreInCellArray = true;
        break
      end
    end
    
    % Increment Dimensions in Acquisition Counter
    for d = nDimensions:-1:1
      if DimensionIndices(d)<Vary.Points(d)
        DimensionIndices(d) = DimensionIndices(d)+1;
        break;
      else
        DimensionIndices(d) = 1;
      end
    end
    
  end
  
  % Reset Acquisition counter to first data point
  DimensionIndices = ones(1,nDimensions);
  
elseif Resonator
  % This is called, if the simulation is a single experiment (no
  % Exp.nPoints and hence no Vary structure)
  for iEvent = 1 : nEvents-1
    if strcmp(Events{iEvent}.type,'pulse')
      if strcmp(Events{iEvent+1}.type,'pulse')
        error('When using a resonator, pulses need to be separated by inter pulse delays to accomodate for ringing from the resonator.')
      end
      % Inter pulse delay of the following event is shortened by the
      % (rounded) duration of ringing
      Ringing = Events{iEvent}.Ringing;
      DecimalRound = 10;
      Ringing = round(Ringing*(10^DecimalRound))/(10^DecimalRound);
      Events{iEvent+1}.t = Events{iEvent+1}.t - Ringing;
    end
    
  end
end

%----------------------------------------------------------------------
% Loop over number of Points = Product of the Points in each
% Dimension
%----------------------------------------------------------------------
for iPoints = 1 : nPoints
  
  EventLengths = InitialEventLengths;
  Sigma = initialSigma;
  
  %--------------------------------------------------------------------
  % Overwrites IQ and t axis for all pulses that are being modified.
  % This needs to be adapted so that only events are overwritten that
  % are actually changed!
  %--------------------------------------------------------------------
  if~isempty(Vary)
    % Replace the inital lengths of the delays with the one that
    % correspond to the current acquisition
    for iDimension = 1 : nDimensions
      if any(any(Vary.IncrementationTable{iDimension}))
        ModifiedEvents = find(Vary.IncrementationTable{iDimension}(:,DimensionIndices(iDimension)));
        for i = 1 : length(ModifiedEvents)
          EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable{iDimension}(ModifiedEvents(i),DimensionIndices(iDimension));
        end
      end
    end
    
    % Load the pulse shapes corresponding to the current acquisition
    % point
    for iPulse = 1 : nPulses
      if ~isempty(Vary.Pulses{iPulse})
        iEvent = PulsePositions(iPulse);
        IQ = LoadfromArray(Vary.Pulses{iPulse}.IQs,DimensionIndices);
        t = LoadfromArray(Vary.Pulses{iPulse}.ts,DimensionIndices);
        % Only replace if the new wave form is different from the
        % previous one
        if ~isequal(Events{iEvent}.IQ,IQ)
          Events{iEvent}.IQ = IQ;
          Events{iEvent}.t = t;
          % Clear propagators
          Events{iEvent}.Propagation = [];
        end
        % Write pulse length to EventLengths vector
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
  
  %--------------------------------------------------------------------
  
  %--------------------------------------------------------------------
  % Check if the Propagators that are needed for this acquisition are
  % available and if not, build and store them
  %--------------------------------------------------------------------
  for iEvent = 1 : nEvents
    currentEvent = Events{iEvent};
    
    if isempty(currentEvent.Propagation) || ~(isfield(currentEvent.Propagation,'Utotal') || isfield(currentEvent.Propagation,'Ltotal'))
      
      switch currentEvent.type
        case 'pulse'
          
          nPhaseCycle = size(currentEvent.PhaseCycle,1);
          %------------------------------------------------------------
          % convert the IQ wave form into a binary form, so that it
          % is possible to use a propagator look up table
          %------------------------------------------------------------
          
          rMax = max(abs(real(currentEvent.IQ(:)))); % max(abs(real....(:)))
          iMax = max(abs(imag(currentEvent.IQ(:))));
          
          MaxWave = max(rMax,iMax);
          
          if isfield(currentEvent,'vertRes')
            % undocumented feature, you can provide your pulses with
            % their own vertical resolution e.g. for a rectangular
            % pulse it would suffice to have a vertRes of 3 (-1,0,1)
            if mod(currentEvent.vertRes,2)
              vertRes = currentEvent.vertRes+1;
            else
              vertRes = currentEvent.vertRes;
            end
          else
            vertRes = 1024;
          end
          
          scale = 2*2*2*pi*MaxWave/vertRes;
          % one factor 2 is required because of linearly polarized
          % irradiation, the other because of the digitization of the
          % wave
          
          tvector(1) = 0;
          tvector(2:length(currentEvent.t)+1) = currentEvent.t+currentEvent.TimeStep;
          
          realArbitrary = real(currentEvent.IQ)/MaxWave;
          realBinary = round(vertRes*(realArbitrary+1)/2);
          
          if currentEvent.ComplexExcitation == 1
            imagArbitrary = imag(currentEvent.IQ)/MaxWave;
            imagBinary = round(vertRes*(imagArbitrary+1)/2);
          end
          
          %------------------------------------------------------------
          % Propagator/Liouvillian Calculation Starts Here
          %------------------------------------------------------------
          % Only build the propagators if they are not available
          if (~currentEvent.Relaxation && ~isfield(currentEvent.Propagation,'Utotal')) ... % checks for availability if simulation is in Hilbert space
              || (currentEvent.Relaxation && ~isfield(currentEvent.Propagation,'Ltotal')) % checks for availability in Liouville space
            
            if currentEvent.Relaxation % Set up for Liouville space simulations
              n = size(Sigma,1);
              Gamma = Relaxation.Gamma;
              equilibriumState = reshape(Relaxation.equilibriumState,n*n,1);
            end
            
            %----------------------------------------------------------
            % Build the Propagator table, only possible for
            % excitation that has no imaginary part - increases speed
            % by precomputing all propagators and later loading them,
            % instead of calculating them at each way point (true as
            % long as length(Wave) > vertical resolution)
            %----------------------------------------------------------
            if currentEvent.ComplexExcitation == 0 && ~currentEvent.Relaxation
              
              UTable = buildPropagators(Ham0,currentEvent.xOp,currentEvent.TimeStep,vertRes,scale);
              
              Events{iEvent}.Propagation.UTable = UTable;
              
            elseif currentEvent.ComplexExcitation == 0 && currentEvent.Relaxation
              
              [LTable, SigmassTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,currentEvent.xOp,currentEvent.TimeStep,vertRes,scale);
              
              Events{iEvent}.Propagation.LTable = LTable;
              Events{iEvent}.Propagation.SigmassTable = SigmassTable;
              
            end
            
            % Set up structures for the actual propagators
            if currentEvent.Detection && ~currentEvent.Relaxation
              Utotal = cell(nPhaseCycle,length(realBinary));
            else
              if ~currentEvent.Relaxation
                Utotal = cell(nPhaseCycle,1);
              else
                Ltotal = cell(nPhaseCycle,length(realBinary));
                SigmaSStotal = Ltotal;
              end
            end
            
            for iPhaseCycle = 1 : nPhaseCycle
              for iWavePoint = 1 : length(realBinary)
                if currentEvent.ComplexExcitation == 0
                  %----------------------------------------------------
                  % Lookup tables are only feasible for non-complex
                  % excitation
                  %----------------------------------------------------
                  if ~currentEvent.Relaxation
                    U = UTable{realBinary(iPhaseCycle,iWavePoint)+1};
                  else
                    L=LTable{realBinary(iPhaseCycle,iWavePoint)+1};
                    SigmaSS=SigmassTable{realBinary(iPhaseCycle,iWavePoint)+1};
                  end
                else
                  %----------------------------------------------------
                  % For complex excitation it is not possible to use a
                  % lookup table
                  %----------------------------------------------------
                  Ham1 = scale/2*(realBinary(iPhaseCycle,iWavePoint)-vertRes/2)*real(currentEvent.xOp)+1i*scale/2*(imagBinary(iPhaseCycle,iWavePoint)-vertRes/2)*imag(currentEvent.xOp);
                  Ham =  Ham0+Ham1;
                  if ~currentEvent.Relaxation
                    U = Propagator(Ham,currentEvent.TimeStep);
                  else
                    [L, SigmaSS] = Liouvillian(Ham,Gamma,equilibriumState,currentEvent.TimeStep);
                  end
                end
                
                if ~currentEvent.Relaxation
                  %------------------------------------------------------
                  % Store Propagators for propation in one step
                  % (additional speed bost if no detection) or stepwise
                  %------------------------------------------------------
                  if currentEvent.Detection
                    Utotal{iPhaseCycle,iWavePoint} = U;
                  else
                    if iWavePoint == 1
                      Utotal{iPhaseCycle,1} = U;
                    else
                      Utotal{iPhaseCycle,1} = U*Utotal{iPhaseCycle,1};
                    end
                  end
                else
                  Ltotal{iPhaseCycle,iWavePoint} = L;
                  SigmaSStotal{iPhaseCycle,iWavePoint} = SigmaSS;
                  
                end
              end
            end
            
            
            %-------------------------------------------------------
            % Write to Event Structure
            %--------------------------------------------------------
            if ~currentEvent.Relaxation
              Events{iEvent}.Propagation.Utotal = Utotal;
            else
              Events{iEvent}.Propagation.Ltotal = Ltotal;
              Events{iEvent}.Propagation.SigmaSStotal = SigmaSStotal;
            end
          end
          
      end
      
    end
  end
  
  %--------------------------------------------------------------------
  % Preparing for the actual propagation step
  %--------------------------------------------------------------------
  % Setting up some initial variables
  
  % bookkeping to keep track of when the first event is detected and
  % what should be stored, needed for proper arrangment of the output
  % later
  ttotal = 0;
  firstDetection = true;
  firstDensityMatrix = true;
  CurrentStateTrajectories = [];
  startTrace = 2;
  
  % Get the data points in all dimensions of the current acquistion,
  % will be used to store the output at the correct position in the
  % multidimensional arrays
  AcquisitionIndex = num2cell(DimensionIndices);
  
  %--------------------------------------------------------------------
  % Loop over event structure
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
          currentSignal = zeros(nDet,size(currentEvent.IQ,2)+1);
        case 'free evolution'
          tvector = 0:currentEvent.TimeStep:currentEvent.t;
          currentSignal = zeros(nDet,length(tvector));
      end
      
      n = size(Sigma,1);
      
      for iDet = 1:length(Det)
        Det{iDet} = reshape(Det{iDet}.',1,n^2);
        normsDet(iDet) = Det{iDet}*Det{iDet}';
      end
      currentSignal(:,1) = Detect(Sigma,Det,normsDet);
      
    else
      switch currentEvent.type
        case 'pulse'
          tvector = currentEvent.t;
        case 'free evolution'
          tvector = [0 currentEvent.t];
          currentEvent.TimeStep = currentEvent.t;
      end
      
      currentSignal=[];
    end
    
    if currentEvent.StateTrajectories
      if currentEvent.Detection
        DensityMatrices=cell(1,length(currentEvent.t));
      else
        DensityMatrices=cell(1,2);
      end
      DensityMatrices{1}=Sigma;
    else
      DensityMatrices = [];
    end
    %------------------------------------------------------------------
    
    switch currentEvent.type
      case 'pulse'
        %--------------------------------------------------------------
        % This creates a new time axis from the one provided inside the
        % event structure. This is necessary, as the first point should
        % represent the state of the system as it was before any
        % irradiation. The pulse therefore starts at t0 + dt;
        %--------------------------------------------------------------
        tvector(1) = 0;
        tvector(2:length(currentEvent.t)+1) = currentEvent.t+currentEvent.TimeStep;
        
        nPhaseCycle = size(currentEvent.PhaseCycle,1);
        
        %--------------------------------------------------------------
        % If phasecycling is requested for this event, a few variables
        % need to be defined
        %--------------------------------------------------------------
        if nPhaseCycle>1
          StateBeforePC = Sigma;
          loopState = zeros(size(Sigma));
          PCnorm = sum(abs(currentEvent.PhaseCycle(:,2)));
        end
        
        %--------------------------------------------------------------
        % Propagation Starts Here, and loops over all steps of the
        % phase cycle
        %--------------------------------------------------------------
        for iPhaseCycle = 1 : nPhaseCycle
          
          if ~currentEvent.Relaxation
            nSteps = size(currentEvent.Propagation.Utotal,2);
          else
            nSteps = size(currentEvent.Propagation.Ltotal,2);
            SigmaVector = reshape(Sigma,n*n,1);
          end
          
          %------------------------------------------------------------
          % Loop over the number of necessary propagation steps (1 if
          % detection during the pulse is not requested
          %------------------------------------------------------------
          for iWavePoint = 1 : nSteps
            if ~currentEvent.Relaxation
              U = currentEvent.Propagation.Utotal{iPhaseCycle,iWavePoint};
              Sigma = U*Sigma*U';
            elseif currentEvent.Relaxation
              L = currentEvent.Propagation.Ltotal{iPhaseCycle,iWavePoint};
              SigmaSS = currentEvent.Propagation.SigmaSStotal{iPhaseCycle,iWavePoint};
              SigmaVector = SigmaSS+L*(SigmaVector-SigmaSS);
              
              Sigma = reshape(SigmaVector,n,n);
            end
            
            %----------------------------------------------------------
            % Store expectation values and density matrix if requested
            %----------------------------------------------------------
            if currentEvent.Detection
              currentSignal(:,iWavePoint+1) = Detect(Sigma,Det,normsDet);
            end
            
            if currentEvent.StateTrajectories
              DensityMatrices{iWavePoint+1} = Sigma;
            end
            %----------------------------------------------------------
          end
          
          %------------------------------------------------------------
          % Combine Results from current phase cycle with previous ones
          %------------------------------------------------------------
          if nPhaseCycle > 1
            PCweight = currentEvent.PhaseCycle(iPhaseCycle,2)/PCnorm;
            loopState = loopState+PCweight*Sigma;
            if currentEvent.Detection
              if iPhaseCycle ~= 1
                PCSignal = PCSignal+PCweight*currentSignal;
              else
                PCSignal = PCweight*currentSignal;
              end
            end
            if iPhaseCycle ~= nPhaseCycle
              Sigma = StateBeforePC;
            end
          end
          %------------------------------------------------------------
          
        end
        
        %--------------------------------------------------------------
        % If phase cycling was requested, the results from combining
        % the individual signals and density matrices are used to
        % overwrite Sigma and currentSignal
        %--------------------------------------------------------------
        if nPhaseCycle > 1
          Sigma = loopState;
          if currentEvent.Detection
            currentSignal = PCSignal;
          end
        end
        %--------------------------------------------------------------
        
      case 'free evolution'
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
            SigmaSS = currentEvent.Propagation.SigmaSStotal;
          else
            Gamma = Relaxation.Gamma;
            equilibriumState = reshape(Relaxation.equilibriumState,n*n,1);
            [L, SigmaSS] = Liouvillian(Ham0,Gamma,equilibriumState,currentEvent.TimeStep);
            Events{iEvent}.Propagation.Ltotal = L;
            Events{iEvent}.Propagation.SigmaSStotal = SigmaSS;
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
            SigmaVector = SigmaSS+L*(SigmaVector-SigmaSS);
            Sigma = reshape(SigmaVector,n,n);
          end
          
          %------------------------------------------------------------
          % Stores expectation values and denisty matrices if requested
          %------------------------------------------------------------
          if currentEvent.Detection
            currentSignal(:,itvector) = Detect(Sigma,Det,normsDet);
          end
          
          if currentEvent.StateTrajectories
            DensityMatrices{itvector} = Sigma;
          end
          %------------------------------------------------------------
        end
    end
    
    %------------------------------------------------------------------
    % This combines the signals and time axes from all detected events
    % and prepares the output
    %------------------------------------------------------------------
    if ~isempty(currentSignal)
      if firstDetection
        Signal = [];
        t = [];
        % store first point of the first signal
        Signal(:,1) = currentSignal(:,1);
        t(1,1) = ttotal;
        
        % now add all the others timepoints
        nSignal = size(currentSignal,2);
        endTrace = startTrace+nSignal-2;
        Signal(:,startTrace:endTrace) = currentSignal(:,2:end);
        t(1,startTrace:endTrace) = tvector(2:end)+ttotal;
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
        t(1,startTrace:endTrace) = tvector(2:end)+ttotal;
        startTrace = endTrace+1;
      end
    end
    
    % If StateTrajectories are required, they are stored in a temporary
    % cell array, which will later be stored to the output
    if currentEvent.StateTrajectories
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
    % Update Total Time, necessary to keep correct timings  if events
    % are not detected
    ttotal = ttotal + tvector(end);
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
        if length(IndirectDimensions) == 1
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
        % double check that the most recent trace has the same length as
        % all the others before
        % in very rare cases it can happen that the total detection time
        % between different acquistions is identical, but the total number
        % of points changes. This is due to the fact, that if the length of
        % an event is changed along one of the indirect dimensions, this
        % change might not scale with the time step (think of it as the
        % least common denominator). An extreme example would be:
        % 1st Acquisition: tau1 = 0.5 us, tau2 = 0.5 us, dt = 0.2 us
        %                  total detection time = 0.5 + 0.5 = 1 us
        %                  total number of data points = 2 + 2 = 4
        % 2nd Acquisition: tau1 = 0.4 us, tau2 = 0.6 us, dt = 0.2 us
        %                  total detection time = 0.4 + 0.6 = 1 us
        %                  total number of data points = 2 + 3 = 5
        % These are usually very minor numerical differences (not
        % necessesarily errors) that come from the time discretization step
        % and are not expected to affect the simulation (especially since
        % the time step is usually very small compared to any time
        % increment along an indirect dimension).
        StoreInCellArray = true;
        
        % Create cell arrays for output
        if length(IndirectDimensions) == 1
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
          [indexID{:}] = ind2sub(cell2mat(IndirectDimensions),jSignal);

          % assign jSignal to cell array, remove singleton dimensions
          NewSignalArray{indexID{:}} = squeeze(SignalArray(indexID{:},:,:));
          
          % If only one detection operator is used, the squeeze above also
          % rotates the signal vector and we have to reverse the rotation.
          if length(Det) == 1
            NewSignalArray{indexID{:}} = NewSignalArray{indexID{:}}.';
          end
          
          % For the time axis, two cases need to be considered as well. If
          % there is only one indirect dimension, the time vector can be
          % taken as is
          if length(IndirectDimensions) == 1
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
  
  % Store the final state at its correct position
  FinalStates(AcquisitionIndex{:},:,:) = Sigma;
  
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


function UTable = buildPropagators(Ham0,xOp,dt,vertRes,scale)

UTable = cell(1,vertRes+1);

for iRes = 0:vertRes
  Ham1 = scale*(iRes-vertRes/2)*real(xOp);
  Ham = Ham0+Ham1;
  U = Propagator(Ham,dt);
  UTable{iRes+1} = U;
end

  function U = Propagator(Ham,dt)
    U = expm(-1i*Ham*dt);
    
    function [LTable, SigmassTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,xOp,dt,vertRes,scale)
      
      SigmassTable = cell(1,vertRes+1);
      LTable = cell(1,vertRes+1);
      
      for ivertRes = 0:vertRes
        Ham1 = scale*(ivertRes-vertRes/2)*real(xOp);
        Ham = Ham0 + Ham1;
        [L, SigmaSS] = Liouvillian(Ham,Gamma,equilibriumState,dt);
        LTable{ivertRes+1} = L;
        SigmassTable{ivertRes+1} = SigmaSS;
      end
      
      function [L, SigmaSS] = Liouvillian(Ham,Gamma,equilibriumState,dt)
        
        n = size(Ham,1);
        
        HamSuOp = kron(eye(n,n),Ham)-kron(Ham.',eye(n,n));
        L = -1i*HamSuOp-Gamma;
        SigmaSS = Gamma*equilibriumState; % steady state solutions for the density matrices
        SigmaSS = -L\SigmaSS;
        L = expm(L*dt); %calculations of the Liouvillians
        
        function detectedSignal = Detect(Sigma,Det,normsDet)
          
          nDet = numel(Det);
          Density = Sigma(:);
          detectedSignal = zeros(nDet,1);
          
          for iDet = 1:nDet
            detectedSignal(iDet) = Det{iDet}*Density/normsDet(iDet);
          end
          
          function LoadedElement = LoadfromArray(Array, ArrayIndex)
            % This function loads an element from an Array (numeric or cell) for a
            % given ArrayIndex (vector)
            nPerDimension = size(Array);
            IndexToLoad = ones(1,length(ArrayIndex));
            if length(ArrayIndex) == 1
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