% thyme  Time domain evolution of density matrix
function [TimeArray, SignalArray, FinalStates, AllDensityMatrices, Events] = thyme(Sigma,Ham0,Det,Events,Relaxation,Vary)

if (nargin==0), help(mfilename); return; end

% if (nargout<1), error('Not enough output arguments.'); end
if (nargout>5), error('Too many output arguments.'); end
if (nargin<4) || (nargin>6), error('Wrong number of input arguments!'); end

%move this this into Events?
% if (nargin<6), IncScheme = 1; end
% if (nargin<7), Mix = {}; end
%
% if any(mod(n,1)~=0) || any(n<=0)
%     error('n, the number of points (4th argument), must be a positive integer.');
% end

method = 'stepwise';

nEvents = length(Events);

% questions:
% - we have to look into the normalization for the Detection operators. +
% or - Detection operators need to divided by two before normalizations are
% computed. what should the normalization be for transition selective
% operators or for detecting populations

switch method
  
  case 'stepwise'
    
    nDet = numel(Det);
    normsDet = zeros(1,nDet);
    Ham0 = Ham0*2*pi;
    
    if ~isempty(Vary)
      nDimensions = numel(Vary.Points);
      idx = ones(1,nDimensions);
      nPoints = prod(Vary.Points);
    else
      nDimensions = 0;
      nPoints = 1;
    end
         
    n = size(Sigma,2);
    FinalStates = zeros(nPoints,n,n);
    AllDensityMatrices = [];
    initialSigma = Sigma;
    
%     timer = 0;
    
%     PulsePositions = [];
%     iPulse = 1;
    isPulse = zeros(1,nEvents);
    iPulse = 1;
    for iEvent = 1 : nEvents
      switch Events{iEvent}.type
        case 'pulse'
          isPulse(iEvent) = 1;
          PulsePositions(iPulse) = iEvent;
          iPulse = iPulse + 1;
      end
    end
    
    %----------------------------------------------------------------------
    % This checks if traces in all dimensions will have the save length, if
    % not, a cell array is used for storage of traces
    %----------------------------------------------------------------------
    StoreInArray = false;
    

    
    if ~isempty(Vary)
      
      InitialEventLengths = zeros(1,nEvents);
      for iEvent = 1 : nEvents
        InitialEventLengths(iEvent) = Events{iEvent}.t(end);
      end
     for iPoint = 1 : nPoints
       
       EventLengths = InitialEventLengths;
       DetectionTime = 0;
       
       for iDimension = 1 : nDimensions
         
         if any(any(Vary.IncrementationTable{iDimension}))
           ModifiedEvents = find(Vary.IncrementationTable{iDimension}(:,idx(iDimension)));
           for i = 1 : length(ModifiedEvents)
             EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable{iDimension}(ModifiedEvents(i),idx(iDimension));
           end
         end
         
       end
        
        nPulses = length(Vary.Pulses);
        
        for iPulse = 1 : nPulses
          if ~isempty(Vary.Pulses{iPulse})
           iEvent = PulsePositions(iPulse);
           t = LoadWaveform(Vary.Pulses{iPulse}.ts,idx);
           EventLengths(iEvent) = t(end);
          end
            
        end
          
        [Sequence, NewEventLengths] = reorder_events(EventLengths,isPulse);
        
        for iSequence = 1 : nEvents
          if Events{Sequence(iSequence)}.Detection
            DetectionTime = DetectionTime + NewEventLengths(iSequence);
          end
        end
        
        if iPoint == 1
          ReferenceDetection = DetectionTime;
        else
          if abs(ReferenceDetection-DetectionTime) > 10^-12
            
           StoreInArray = true;
           break
           
          end
        end
        
        
        for d = nDimensions:-1:1
          if idx(d)<Vary.Points(d)
            idx(d) = idx(d)+1;
            break;
          else
            idx(d) = 1;
          end
        end
        
     end
     
     idx = ones(1,nDimensions);
     
    end
        
    if StoreInArray
      SignalArray = cell(1,nPoints);
      TimeArray = SignalArray;
    end
    %----------------------------------------------------------------------
    
    
      
%     
%     ModifiedDelays = [];
%     SwappedPulses = [];
    %----------------------------------------------------------------------
    % Loop over number of Points = Product of the Points in each
    % Dimension
    %----------------------------------------------------------------------
    for iPoints = 1 : nPoints
      tic
      
      EventLengths = InitialEventLengths;
            
%       if ~isempty(ModifiedDelays)
%         for iEvent = ModifiedDelays 
%           Events{iEvent} = OrigEvents{iEvent};
%         end
%         
%         PulseSwap = Events{SwappedPulses(1)};
%         Events{SwappedPulses(1)} = Events{SwappedPulses(2)};
%         Events{SwappedPulses(2)} = PulseSwap;
%       end

      Sigma = initialSigma;
      
      %--------------------------------------------------------------------
      % Overwrites IQ and t axis for all events that are being modified
      % This needs to be adapted so that only events are overwritten that
      % are actually changed!
      %--------------------------------------------------------------------
      
      if~isempty(Vary)
                
        for iDimension = 1 : nDimensions
          
          if any(any(Vary.IncrementationTable{iDimension}))
            ModifiedEvents = find(Vary.IncrementationTable{iDimension}(:,idx(iDimension)));
            for i = 1 : length(ModifiedEvents)
              EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable{iDimension}(ModifiedEvents(i),idx(iDimension));
            end
          end
          
        end
        
        nPulses = length(Vary.Pulses);
        
        for iPulse = 1 : nPulses
          if ~isempty(Vary.Pulses{iPulse})
           iEvent = PulsePositions(iPulse);
           Events{iEvent}.IQ = LoadWaveform(Vary.Pulses{iPulse}.IQs,idx);
           Events{iEvent}.t = LoadWaveform(Vary.Pulses{iPulse}.ts,idx);
           EventLengths(iEvent) = Events{iEvent}.t(end);
           Events{iEvent}.Propagation = [];
          end
            
        end
          
        [Sequence, NewEventLengths] = reorder_events(EventLengths,isPulse);
        
        for iEvent = 1 : nEvents
          NewPosition = find(Sequence == iEvent);
          if strcmp(Events{iEvent}.type,'free evolution') && NewEventLengths(NewPosition) ~= Events{iEvent}.t
            Events{iEvent}.t = NewEventLengths(NewPosition);
            Events{iEvent}.Propagation = [];
          end
        end
          
      else
        Sequence = 1 : nEvents;        
      end

      %--------------------------------------------------------------------
           
      %--------------------------------------------------------------------
      % Check if the Propagators that are needed for this acquisition are
      % available and if not, build and store them
      %--------------------------------------------------------------------
      for iEvent = 1 : nEvents
        currentEvent = Events{iEvent};
        
        if isempty(currentEvent.Propagation) || ~((isfield(currentEvent.Propagation,'Utotal') || isfield(currentEvent.Propagation,'Ltotal')))
                    
          switch currentEvent.type
            case 'pulse'
                          
              nPhaseCycle = size(currentEvent.PhaseCycle,1);
              %------------------------------------------------------------
              % convert the IQ wave form into a binary form, so that it
              % is possible to use a propagator look up table
              %------------------------------------------------------------
              
              rMax = max(abs(real(currentEvent.IQ(:)))); % max(abs(real....(:)))
              iMax = max(abs(imag(currentEvent.IQ(:))));
              
              if rMax > iMax
                MaxWave = rMax;
              else
                MaxWave = iMax;
              end
              
              vertRes = 1024;
              
              scale = 2*2*2*pi*MaxWave/vertRes;
              % one factor 2 is required because of linearly polarized 
              % irradiation, the other because of the digitization of the 
              % wave
              
              tvector(1) = 0;
              tvector(2:length(currentEvent.t)+1) = currentEvent.t+currentEvent.TimeStep;
              
              realArbitrary = real(currentEvent.IQ)/MaxWave;
              realBinary = floor(vertRes*(realArbitrary+1)/2);
              realBinary(realBinary == vertRes) = vertRes-1;
              
              if currentEvent.ComplexExcitation == 1
                imagArbitrary = imag(currentEvent.IQ)/MaxWave;
                imagBinary = floor(vertRes*(imagArbitrary+1)/2);
                imagBinary(imagBinary==vertRes) = vertRes-1;
              end
              
              %------------------------------------------------------------
              % Propagator/Liouvillian Calculation Starts Here
              %------------------------------------------------------------
              if ~currentEvent.Relaxation
                %----------------------------------------------------------
                % Load or build the Propagator table, only possible for
                % excitation that has no imaginary part
                %----------------------------------------------------------
                if currentEvent.ComplexExcitation == 0

                  if ~isempty(currentEvent.Propagation) && isfield(currentEvent.Propagation,'UTable')
                    UTable = currentEvent.Propagation.UTable;
                  else
                    UTable = buildPropagators(Ham0,currentEvent.xOp,currentEvent.TimeStep,vertRes,scale);
                    Events{iEvent}.Propagation.UTable = UTable;
                  end
                end
                %----------------------------------------------------------
                
                if currentEvent.Detection
                  Utotal = cell(nPhaseCycle,length(realBinary));
                else
                  Utotal = cell(nPhaseCycle,1);
                end
                
                %----------------------------------------------------------
                % Loop over the Phasecycles
                for iPhaseCycle = 1 : nPhaseCycle
                  % Propagation for one waveform
                  for iWavePoint = 1 : length(realBinary)
                    if currentEvent.ComplexExcitation == 0
                      % Load propagators if Complex Excitation is off
                      U = UTable{realBinary(iPhaseCycle,iWavePoint)+1};
                    else % For active Complex Excitation Propagators need to be recalculated
                      Ham1 = scale/2*(realBinary(iPhaseCycle,iWavePoint)-vertRes/2)*real(currentEvent.xOp)+1i*scale/2*(imagBinary(iPhaseCycle,iWavePoint)-vertRes/2)*imag(currentEvent.xOp);
                      Ham =  Ham0+Ham1;
                      U = Propagator(Ham,currentEvent.TimeStep);
                    end
                    
                    %------------------------------------------------------
                    %  Store Propagators for propation in one step or
                    %  stepwise
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
                    %------------------------------------------------------
                  end
                end
                
                %--------------------------------------------------------
                % Write to Event Structure
                %--------------------------------------------------------
                Events{iEvent}.Propagation.Utotal = Utotal;
                              
              elseif currentEvent.Relaxation
                
                n = size(Sigma,1);
                Gamma = Relaxation.Gamma;
                equilibriumState = reshape(Relaxation.equilibriumState,n*n,1);
                
                %----------------------------------------------------------
                % Check wether a Liouvillian table is available and build
                % one if not
                %----------------------------------------------------------
                if ~isfield(currentEvent.Propagation,'Ltotal')
                  if currentEvent.ComplexExcitation == 0
                    if ~isempty(currentEvent.Propagation) && isfield(currentEvent.Propagation,'LTable')
                      LTable = currentEvent.Propagation.LTable;
                      SigmassTable = currentEvent.Propagation.SigmassTable;
                    else
                      
                      [LTable, SigmassTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,currentEvent.xOp,currentEvent.TimeStep,vertRes,scale);
                      
                      Events{iEvent}.Propagation.LTable = LTable;
                      Events{iEvent}.Propagation.SigmassTable = SigmassTable;
                    end
                    
                  end
                  
                  if currentEvent.Detection
                    Ltotal = cell(nPhaseCycle,length(realBinary));
                    SigmaSStotal = Ltotal;
                  end
                  
                  %--------------------------------------------------------
                  % Calculating Propagators for all steps of phase cycling
                  %--------------------------------------------------------
                  
                  for iPhaseCycle = 1 : nPhaseCycle
                    for iWavePoint = 1 : length(realBinary)
                      %----------------------------------------------------
                      % For complex excitation it is not possible to use a
                      % lookup table
                      %----------------------------------------------------
                      if currentEvent.ComplexExcitation == 0
                        L=LTable{realBinary(iPhaseCycle,iWavePoint)+1};
                        SigmaSS=SigmassTable{realBinary(iPhaseCycle,iWavePoint)+1};
                      else
                        Ham1 = scale/2*(realBinary(iPhaseCycle,iWavePoint)-vertRes/2)*real(currentEvent.xOp)+1i*scale/2*(imagBinary(iPhaseCycle,iWavePoint)-vertRes/2)*imag(currentEvent.xOp);
                        Ham = Ham0+Ham1;
                        [L, SigmaSS] = Liouvillian(Ham,Gamma,equilibriumState,currentEvent.TimeStep);
                      end
                      
                      Ltotal{iPhaseCycle,iWavePoint} = L;
                      SigmaSStotal{iPhaseCycle,iWavePoint} = SigmaSS;
                    end
                  end
                  %--------------------------------------------------------
                  % Write to Event Structure
                  %--------------------------------------------------------
                  Events{iEvent}.Propagation.Ltotal = Ltotal;
                  Events{iEvent}.Propagation.SigmaSStotal = SigmaSStotal;
                end
              end
              
          end
          
        end
      end
      
      
      %--------------------------------------------------------------------
      % Setting up some initial variables
      ttotal = 0;
      firstDetection = true;
      firstDensityMatrix = true;
      startTrace = 2;
      
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
            currentSignal(:,1) = Detect(Sigma,Det,normsDet);
          end
          
        else
          switch currentEvent.type
            case 'pulse'
              tvector = currentEvent.t;
            case 'free evolution'
              tvector = [0 currentEvent.t];
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
                DensityMatrices{itvector}=Sigma;
              end
              %------------------------------------------------------------
            end            
        end
        
        %------------------------------------------------------------------
        % This combines the signals and time axes from all detected events
        % and prepares the output
        %------------------------------------------------------------------        
        if~isempty(currentSignal)
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
            % adding other signals, this confusing index is necessary in
            % order to avoid double counting last point of a signal and the
            % first point of the succiding signal
            nSignal = size(currentSignal,2);
            endTrace = startTrace+nSignal-2;
            Signal(:,startTrace:endTrace) = currentSignal(:,2:end);
            t(1,startTrace:endTrace) = tvector(2:end)+ttotal;
            startTrace = endTrace+1;
          end
        end
        
        if currentEvent.StateTrajectories
          if firstDensityMatrix
            AllDensityMatrices{iPoints} = DensityMatrices;
            firstDensityMatrix = false;
          else
            iStart = length(AllDensityMatrices{iPoints});
            nElements = length(DensityMatrices);
            for iDensity = 2 : nElements
              AllDensityMatrices{iPoints}{iStart+iDensity-1} = DensityMatrices{iDensity};
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
      % array (cell array if eventlengths changed) that contains all the
      % time traces according to the Vary structure
      %--------------------------------------------------------------------
      if ~firstDetection
        % move all the detected signals into large array here!
        if StoreInArray
          SignalArray{iPoints} = Signal;
          TimeArray{iPoints} = t;
        else
          if iPoints == 1
            SignalSize = size(Signal);
            SignalArray = zeros(nPoints,SignalSize(1),SignalSize(2));
            SignalArray(iPoints,:,:) = Signal;
            TimeArray = zeros(nPoints,SignalSize(2));
            TimeArray(1,:) = t;
          else
            SignalArray(iPoints,:,:) = Signal;
            TimeArray(iPoints,:) = t;
          end
        end
      end
      toc
      
      FinalStates(iPoints,:,:) = Sigma;
      %--------------------------------------------------------------------
      % Incremeant the index for the Vary structure by 1
      %--------------------------------------------------------------------
      
      if ~isempty(Vary)
        for d = nDimensions:-1:1
          if idx(d)<Vary.Points(d)
            idx(d) = idx(d)+1;
            break;
          else
            idx(d) = 1;
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
    %----------------------------------------------------------------------   
    
  case 'incrementation scheme'
    
    
    
    % IncScheme check
    %------------------------------------------------------------
    if (length(IncScheme)>1) && (nargin<7),
      error('The requested IncScheme requires mixing propagators, but none are provided!');
    end
    if any((abs(IncScheme)~=1) & (abs(IncScheme)~=2))
      error('IncScheme can contain only 1, -1, 2, and -2.');
    end
    
    nEvolutionPeriods = numel(IncScheme);
    nDimensions = max(abs(IncScheme));
    
    % Parameter parsing
    %------------------------------------------------------------
    if ~iscell(Det)
      Det = {Det};
    end
    nDetectors = numel(Det);
    
    if ~iscell(Mix)
      Mix = {Mix};
    end
    nMixingBlocks = numel(Mix);
    
    if (nMixingBlocks~=nEvolutionPeriods-1),
      error('Number of mixing propagators not correct! %d are needed.',nEvolutionPeriods-1);
    end
    N = size(Sigma,1);
    
    if (nDimensions==1)
      for iDet = 1:nDetectors
        Signal{iDet} = zeros(n,1);
      end
      if iscell(Ham0), Ham0 = Ham0{1}; end
    else
      if numel(dt)==1, dt = [dt dt]; end
      if numel(n)==1, n = [n n]; end
      for iDet = 1:nDetectors
        Signal{iDet} = zeros(n);
      end
    end
    if nDetectors==1
      Signal = Signal{1};
    end
    
    % Transform all operators to Hamiltonian eigenbasis (eigenbases)
    %---------------------------------------------------------------
    if ~iscell(Ham0)
      
      if nnz(Ham0)==nnz(diag(Ham0)) % Check if Hamiltonian is already diagonal
        E = diag(Ham0);
        Density = Sigma;
        Detector = Det;
      else
        % Diagonalize Hamiltonian
        [Vecs,E] = eig(Ham0); % MHz, E doesn't have to be sorted
        E = real(diag(E));
        % Transform all other matrices to Hamiltonian eigenbasis
        Density = Vecs'*Sigma*Vecs;
        for iMix = 1:nMixingBlocks
          Mix{iMix} = Vecs'*Mix{iMix}*Vecs;
        end
        for iDet = 1:nDetectors
          Detector{iDet} = Vecs'*Det{iDet}*Vecs;
        end
      end
      
      % Define free evolution propagators
      if (nDimensions==1)
        diagU = exp(-2i*pi*dt*E);
      else
        diagUX = exp(-2i*pi*dt(1)*E);
        diagUY = exp(-2i*pi*dt(2)*E);
      end
      
    else
      
      % Check if Hamiltonians are already diagonal
      if (nnz(Ham0{1})==nnz(diag(Ham0{1})) && nnz(Ham0{1})==nnz(diag(Ham0{1})))
        Ex = Ham0{1};
        Ey = Ham0{2};
        Density = Sigma;
        Detector = Det;
      else
        
        % Diagonalize Hamiltonians
        [Vecs{1},Ex] = eig(Ham0{1});
        [Vecs{2},Ey] = eig(Ham0{2});
        
        % Transform all other matrices to Hamiltonian eigenbasis
        d = abs(IncScheme);
        Density = Vecs{d(1)}'*Sigma*Vecs{d(1)};
        for iMix = 1:nMixingBlocks
          Mix{iMix} = Vecs{d(iMix+1)}'*Mix{iMix}*Vecs{d(iMix)};
        end
        for iDet = 1:nDetectors
          Detector{iDet} = Vecs{d(end)}'*Det{iDet}*Vecs{d(end)};
        end
        
      end
      
      % Define free evolution propagators
      diagUX = exp((-2i*pi*dt(1))*real(diag(Ex)));
      diagUY = exp((-2i*pi*dt(2))*real(diag(Ey)));
      
    end
    
    % Time-domain evolution, IncScheme switchyard
    %------------------------------------------------------------
    % The following implementations in the propagator eigenframes
    % (giving diagonal propagators) make use of the following simplifications
    % of the matrix multiplications associated with the propagations:
    %
    %   U*Density*U'   = (diagU*diagU').*Density
    %   U*Propagator*U = (diagU*diagU.').*Propagator
    %   (U^-1)*Propagator*U = U'*Propagator*U = (conj(diagU)*diagU.').*Propagator
    %   U*Propagator*(U^-1) = U*Propagator*U' = (diagU*diagU').*Propagator
    %
    % where U are diagonal matrices and diagU are vectors of eigenvalues.
    
    % Pre-reshape for trace calculation
    for iDet = 1:nDetectors
      Detector{iDet} = reshape(Detector{iDet}.',1,N^2);
    end
    if nDetectors==1
      Detector = Detector{1};
    end
    
    if isequal(IncScheme,1) % IncScheme [1]
      FinalDensity = Density(:);
      U_ = diagU*diagU';
      U_ = U_(:);
      for ix = 1:n
        % Compute trace(Detector*FinalDensity)
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity;
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity;
          end
        end
        FinalDensity = U_.*FinalDensity; % equivalent to U*FinalDensity*U'
      end
      
    elseif isequal(IncScheme,[1 1]) % IncScheme [1 1]
      UU_ = diagU*diagU.';
      % It is not necessary to evolve the initial density matrix.
      % Only the mixing propagator needs to be evolved.
      Mix1 = Mix{1};
      for ix = 1:n
        % Compute density right before Detection
        FinalDensity = Mix1*Density*Mix1';
        % Compute trace(Detector*FinalDensity)
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1 = UU_.*Mix1; % equivalent to U*Mix1*U
      end
      
    elseif isequal(IncScheme,[1 -1]) % IncScheme [1 -1]
      %   % Pre-propagate mixing propagator to end of second period (= start of
      %   % experiment)
      %   MixX = diag(diagU.^n)*Mix{1};
      MixX =Mix{1};
      UtU_ = conj(diagU)*diagU.';
      for ix = 1:n
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        MixX = UtU_.*MixX; % equivalent to U^-1*MixX*U
      end
      
    elseif isequal(IncScheme,[1 2]) % IncScheme [1 2]
      UX_ = diagUX*diagUX';
      UY_ = diagUY*diagUY';
      UY_ = reshape(UY_,N^2,1);
      Mix1 = Mix{1};
      for ix = 1:n(1)
        FinalDensity = reshape(Mix1*Density*Mix1',N^2,1);
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity;
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity;
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*Density*UY';
        end
        Density = UX_.*Density; % equivalent to UX*Density*UX';
      end
      
    elseif isequal(IncScheme,[1 1 2]) % IncScheme [1 1 2]
      UUX_ = diagUX*diagUX.';
      UY_ = diagUY*diagUY';
      Mix1 = Mix{1};
      Mix2 = Mix{2};
      for ix = 1:n(1)
        M = Mix2*Mix1;
        FinalDensity = M*Density*M';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        Mix1 = UUX_.*Mix1; % equivalent to UX*Mix1*UX
      end
      
    elseif isequal(IncScheme,[1 -1 2]) % IncScheme [1 -1 2]
      UtUX_ = conj(diagUX)*diagUX.';
      UY_ = diagUY*diagUY';
      %   % Pre-propagate mixing propagator to end of second period (= start of
      %   % experiment)
      %   MixX = diag(diagUX.^n(1))*Mix{1};
      MixX = Mix{1};
      Mix2 = Mix{2};
      for ix = 1:n(1)
        M = Mix2*MixX;
        FinalDensity = M*Density*M';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        MixX = UtUX_.*MixX; % equivalent to U^-1*MixX*U
      end
      
    elseif isequal(IncScheme,[1 2 1]) % IncScheme [1 2 1]
      Mix1 = Mix{1};
      Mix2 = Mix{2};
      UUX_ = diagUX*diagUX.';
      UY = diag(diagUY);
      for iy = 1:n(2)
        MixY = Mix2*Mix1;
        MixYadj = MixY';
        for ix = 1:n(1)
          FinalDensity = MixY*Density*MixYadj;
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
        end
        Mix1 = UY*Mix1;
      end
      
    elseif isequal(IncScheme,[1 2 2 1]) % IncScheme [1 2 2 1]
      Mix1 = Mix{1};
      Mix2 = Mix{2};
      Mix3 = Mix{3};
      UUX_ = diagUX*diagUX.';
      UUY_ = diagUY*diagUY.';
      for iy = 1:n(2)
        MixY = Mix3*Mix2*Mix1;
        MixYadj = MixY';
        for ix = 1:n(1)
          FinalDensity = MixY*Density*MixYadj;
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
        end
        Mix2 = UUY_.*Mix2; % equivalent to UY*Mix2*UY
      end
      
    elseif isequal(IncScheme,[1 2 -2 1]) % IncScheme [1 2 -2 1]
      Mix1 = Mix{1};
      %   Mix2 = diag(diagUY.^n(2))*Mix{2}; % pre-propagate to endpoint of third delay
      Mix2 = Mix{2};
      Mix3 = Mix{3};
      UUX_ = diagUX*diagUX.';
      UtUY_ = conj(diagUY)*diagUY.';
      for iy = 1:n(2)
        MixY = Mix3*Mix2*Mix1;
        MixYadj = MixY';
        for ix = 1:n(1)
          FinalDensity = MixY*Density*MixYadj;
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          MixY = UUX_.*MixY; % equivalent to UX*MixY*UX
        end
        Mix2 = UtUY_.*Mix2; % equivalent to UY'*Mix2*UY
      end
      
    elseif isequal(IncScheme,[1 -1 1 -1]) % IncScheme [1 -1 1 -1]
      %   Mix1X = diag(diagU.^n)*Mix{1}; % pre-propagate to endpoint of second delay
      Mix1X = Mix{1};
      Mix2 = Mix{2};
      %   Mix3X = diag(diagU.^n)*Mix{3}; % pre-propagate to endpoint of fourth delay
      Mix3X = Mix{3};
      UtU_ = conj(diagU)*diagU.'; % propagator for Mix1 and Mix3 (add before, remove after)
      for ix = 1:n
        MixX = Mix3X*Mix2*Mix1X;
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1X = UtU_.*Mix1X; % equivalent to U'*Mix1X*U
        Mix3X = UtU_.*Mix3X; % equivalent to U'*Mix3X*U
      end
      
    elseif isequal(IncScheme,[1 1 -1 -1]) % IncScheme [1 1 -1 -1]
      Mix1 = Mix{1};
      %   Mix2X = diag(diagU.^n)*Mix{2}; % pre-propagate to endpoint of third delay
      %   Mix3X = diag(diagU.^n)*Mix{3}; % pre-propagate to endpoint of fourth delay
      Mix2X = Mix{2};
      Mix3X = Mix{3};
      UU1_ = diagU*diagU.'; % propagator for Mix1 (add before and after)
      UU3_ = conj(diagU*diagU.'); % propagator for Mix3 (remove before and after)
      for ix = 1:n
        MixX = Mix3X*Mix2X*Mix1;
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1 = UU1_.*Mix1; % equivalent to U*Mix1*U
        Mix3X = UU3_.*Mix3X; % equivalent to U'*Mix3X*U'
      end
      
    elseif isequal(IncScheme,[1 1 -1 -1 2]) % IncScheme [1 1 -1 -1 2]
      Mix1 = Mix{1};
      %   Mix2X = diag(diagUX.^n(1))*Mix{2}; % pre-propagate to endpoint of third delay
      %   Mix3X = diag(diagUX.^n(1))*Mix{3}; % pre-propagate to endpoint of fourth delay
      Mix2X = Mix{2};
      Mix3X = Mix{3};
      Mix4 = Mix{4};
      UU1_ = diagUX*diagUX.'; % propagator for Mix1 (add before and after)
      UU3_ = conj(diagUX*diagUX.'); % propagator for Mix3 (remove before and after)
      UY_ = diagUY*diagUY';
      for ix = 1:n(1)
        MixX = Mix4*Mix3X*Mix2X*Mix1;
        FinalDensity = MixX*Density*MixX';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        Mix1 = UU1_.*Mix1; % equivalent to U*Mix1*U
        Mix3X = UU3_.*Mix3X; % equivalent to U'*Mix3X*U'
      end
      
    elseif isequal(IncScheme,[1 -1 -1 1]) % IncScheme [1 -1 -1 1]
      %   Mix1X = diag(diagU.^n)*Mix{1}; % pre-propagate to endpoint of second delay
      Mix1X = Mix{1};
      Mix2 = Mix{2};
      %   Mix3X = Mix{3}*diag(diagU.^n); % forward-propagate to endpoint of fourth delay
      Mix3X = Mix{3};
      UtU1_ = conj(diagU)*diagU.'; % propagator for Mix1 (add before, remove after)
      UtU3_ = diagU*diagU'; % propagator for Mix3 (remove before, add after)
      for ix = 1:n
        MixX = Mix3X*Mix2*Mix1X;
        FinalDensity = MixX*Density*MixX';
        if nDetectors==1
          Signal(ix) = Detector*FinalDensity(:);
        else
          for iDet = 1:nDetectors
            Signal{iDet}(ix) = Detector{iDet}*FinalDensity(:);
          end
        end
        Mix1X = UtU1_.*Mix1X; % equivalent to U'*Mix1X*U
        Mix3X = UtU3_.*Mix3X; % equivalent to U*Mix3X*U'
      end
      
    elseif isequal(IncScheme,[1 -1 -1 1 2]) % IncScheme [1 -1 -1 1 2]
      %   Mix1X = diag(diagUX.^n(1))*Mix{1}; % pre-propagate to endpoint of second delay
      Mix1X = Mix{1};
      Mix2 = Mix{2};
      %   Mix3X = Mix{3}*diag(diagUX.^n(1)); % forward-propagate to endpoint of fourth delay
      Mix3X = Mix{3};
      Mix4 = Mix{4};
      UtU1_ = conj(diagUX)*diagUX.'; % propagator for Mix1 (add before, remove after)
      UtU3_ = diagUX*diagUX'; % propagator for Mix3 (remove before, add after)
      UY_ = diagUY*diagUY';
      for ix = 1:n(1)
        MixX = Mix4*Mix3X*Mix2*Mix1X;
        FinalDensity = MixX*Density*MixX';
        for iy = 1:n(2)
          if nDetectors==1
            Signal(ix,iy) = Detector*FinalDensity(:);
          else
            for iDet = 1:nDetectors
              Signal{iDet}(ix,iy) = Detector{iDet}*FinalDensity(:);
            end
          end
          FinalDensity = UY_.*FinalDensity; % equivalent to UY*FinalDensity*UY'
        end
        Mix1X = UtU1_.*Mix1X; % equivalent to U'*Mix1X*U
        Mix3X = UtU3_.*Mix3X; % equivalent to U*Mix3X*U'
      end
      
    else
      error('Unsupported incrementation scheme!');
    end
    
    return
end


function UTable = buildPropagators(Ham0,xOp,dt,vertRes,scale)

UTable = cell(1,vertRes);

for iRes = 0:vertRes-1
  Ham1 = scale*(iRes-vertRes/2)*real(xOp);
  Ham = Ham0+Ham1;
  U = Propagator(Ham,dt);
  UTable{iRes+1} = U;
end

function U = Propagator(Ham,dt)
U = expm_fastc(-1i*Ham*dt);

function [LTable, SigmassTable] = buildLiouvillians(Ham0,Gamma,equilibriumState,xOp,dt,vertRes,scale)

SigmassTable = cell(1,vertRes);
LTable = cell(1,vertRes);

for ivertRes = 0:vertRes-1
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
L = expm_fastc(L*dt); %calculations of the Liouvillians

function detectedSignal = Detect(Sigma,Det,normsDet)

nDet = numel(Det);
Density = Sigma(:);
detectedSignal = zeros(nDet,1);

for iDet = 1:nDet
  detectedSignal(iDet) = Det{iDet}*Density/normsDet(iDet);
end

function LoadedElement = LoadWaveform(Array, ArrayIndex)
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
  
  LoadedElement = Array{IndexToLoad{:}};