function [Events, Vary, Opt] = sequencer(Exp,Opt)
% This function creates the Event and Vary Structures

% -------------------------------------------------------------------------
% Pre-Processing
% -------------------------------------------------------------------------
Vary = [];

% Check if Timestep is sufficient
MaxFreq = max(abs(Exp.Frequency(:)));
Nyquist = 2*MaxFreq;

if Exp.TimeStep > 1/Nyquist
  warning('Your Time Step (Exp.TimeStep) appears to not fullfill the Nyquist criterium for the pulses you provided.')
end

% Create an empty cell array for all the events
Events = cell(1,length(Exp.t));
nEvents = length(Exp.t);

% Variables for bookkeeping of pulses and free evolution events
iPulse = 0;
iDelay = 0;

% A vector to quickly identify pulses, required for the reordering if
% pulses cross during the sequence
isPulse = zeros(1,length(Events));

% Setting up data structures for the pulses and events
for iEvent = 1 : nEvents
  if iEvent > length(Exp.Pulses) || ~isstruct(Exp.Pulses{iEvent})
    iDelay = iDelay + 1;
    DelayIndices(iDelay) = iEvent;
  else
    iPulse = iPulse + 1;
    PulseIndices(iPulse) = iEvent;
    isPulse(iEvent) = 1;
  end
end
 
nPulses = length(PulseIndices);
Pulses = cell(1,nPulses);
iPulse = 1;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create the Eventstructure
% -------------------------------------------------------------------------
for iEvent = 1 : length(Exp.t)
  if length(Exp.Pulses) >= iEvent && isstruct(Exp.Pulses{iEvent})
    % ---------------------------------------------------------------------
    % Pulse Specific Fields 
    % ---------------------------------------------------------------------
    Pulse = Exp.Pulses{iEvent};
    
    Pulse.tp = Exp.t(iEvent);
    
    % Gets the frequency band from Exp.Frequency. If only one frequency
    % band is provided it is used for all pulses.
    if size(Exp.Frequency,1) == 1
      Pulse.Frequency = Exp.Frequency;
    elseif size(Exp.Frequency,1) < iPulse
      error('The Frequency Band for Pulse No. %d is missing.',iPulse)
    else
      Pulse.Frequency = Exp.Frequency(iPulse,:);
    end
       
    % Gets the flip angle 
    if length(Exp.Flip) < iPulse
      error('No Flipangle for Pulse No. %d provided.',iPulse)
    else
      Pulse.Flip = Exp.Flip(iPulse);
    end
    
    % Gets the phase for the pulse, if none is provided, the phase is
    % assumed to be 0
    if isfield(Exp,'Phase') && length(Exp.Phase) >= iPulse
      Pulse.Phase = Exp.Phase(iPulse);
    else
      Pulse.Phase = 0;
    end
    
    % Gets the PhaseCycle for the current Pulse, if none is provided, phase
    % cycling is switched off for this event
    if isfield(Exp,'PhaseCycle') &&  iPulse <= length(Exp.PhaseCycle) && ~isempty(Exp.PhaseCycle{iPulse})
      Pulse.PhaseCycle = Exp.PhaseCycle{iPulse};
    else
      Pulse.PhaseCycle = 0;
    end
       
    % Get the time step size if available
    Pulse.TimeStep = Exp.TimeStep;
       
    % Specify Type in Event structure
    Events{iEvent}.type = 'pulse';
    
    % Loop over the function that creates the PulseShape, as many times at
    % are necessary to calculate all wave forms for the phase cycling
    for iPCstep = 1 : size(Pulse.PhaseCycle,1)
      Pulse.Phase = Pulse.Phase+Pulse.PhaseCycle(iPCstep,1);
      [t,IQ] = pulse(Pulse);
      Events{iEvent}.IQ(iPCstep,:) = IQ;
    end
    
    % Store the time axis of the pulse in the Event structure    
    Events{iEvent}.t = t;
    
    % Store the PhaseCycle in the Event structure
    Events{iEvent}.PhaseCycle = Pulse.PhaseCycle;
       
    % Checks if ComplexExcitation is requested for this Pulse, if not
    % specified Complex Excitation is switched off by default - the
    % excitation operator is being built outside of sequencer
    if ~isfield(Pulse,'ComplexExcitation')
      Events{iEvent}.ComplexExcitation = false;
    else
      Events{iEvent}.ComplexExcitation = true;
    end
    
    % Temporarily store pulse paramaters to avoid reassigning them for creating the
    % vary table
    Pulse.EventIndex = iEvent;
    Pulses{iPulse} = Pulse;
    iPulse = iPulse + 1;
  else
    % ---------------------------------------------------------------------
    % Delay/Free Evolution Specific Fields
    % ---------------------------------------------------------------------
    Events{iEvent}.type = 'free evolution';
    Events{iEvent}.t = Exp.t(iEvent);
  end
  
  % -----------------------------------------------------------------------
  % General Fields
  % -----------------------------------------------------------------------
  % The following fields need to be defined for both, pulses and free
  % evolution events
  
  % Check if Relaxation is requested for Events, by default, Relaxation is
  % switched off
  if ~isfield(Opt,'Relaxation')
    Events{iEvent}.Relaxation = false;
  else
    if length(Opt.Relaxation) == 1
      Events{iEvent}.Relaxation = Opt.Relaxation;
    elseif iEvent > length(Opt.Relaxation)
      Events{iEvent}.Relaxation = false;
    else
      Events{iEvent}.Relaxation = Opt.Relaxation(iEvent);
    end
  end
  
  % Check if detection is provided, if no detection is requested, the last
  % event is detected by default
  if ~isfield(Opt,'DetectedEvents')
    if iEvent == length(Exp.t)
      Events{iEvent}.Detection = true;
    else
      Events{iEvent}.Detection = false;
    end
  else
    if length(Opt.DetectedEvents) == 1
      Events{iEvent}.Detection = Opt.DetectedEvents;
    elseif iEvent > length(Opt.DetectedEvents)
      error('You did not specify detection after Event %d.',iEvent);
    else 
      Events{iEvent}.Detection = Opt.DetectedEvents(iEvent);
    end   
  end
  
  % Check if Density Matrices are to be stored, if not specified, Density
  % Matrices are not stored
  if ~isfield(Opt,'StateTrajectories')
    Events{iEvent}.StateTrajectories = false;
  else
    if length(Opt.StateTrajectories) == 1
      Events{iEvent}.StateTrajectories = Opt.StateTrajectories;
    elseif iEvent > length(Opt.StateTrajectories)
      Events{iEvent}.StateTrajectories = false;
    else
      Events{iEvent}.StateTrajectories = Opt.StateTrajectories(iEvent);
    end
  end
  
  % Store an empty propagation structure, will be overwritten by thyme
  Events{iEvent}.Propagation = [];
  
  % Store the time step, which will be needed in thyme to calculate time
  % axis and propagators
  Events{iEvent}.TimeStep = Exp.TimeStep;
    
end

% -------------------------------------------------------------------------
% Creates the Vary structure
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switching incrementation scheme to off for now, this will be an option
% later on, and will decide on how the incremenation tables are stored. For
% the incrementationscheme only linear increments can be used, and the data
% structure can therefore be reduced
IncrementationScheme = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(Exp,'nPoints')
  nDimensions = length(Exp.nPoints);
  Vary.Points = Exp.nPoints;
  
  % Layout the Incrementation data structure
  if IncrementationScheme 
    Vary.IncrementationTable = zeros(nEvents,nDimensions);
  else
    Vary.IncrementationTable = cell(1,nDimensions);
  end
  
  % This cell array will carry all the modifications that are being made to
  % a pulse in any dimension - each element corresponds to a pulse and will
  % contain the dimensions, fields and values that are to be changed 
  PulseModifications = cell(1,nPulses);
  
  % -----------------------------------------------------------------------
  % Loop over all the provided dimensions and process the input - for
  % changing delays or pulse positions, incrementation tables can be
  % created for each dimension. If pulse parameters are changed, the
  % changed pulse parameters are stored in PulseModifications. After all
  % dimensions have been checked, the values stored in Pulse
  % -----------------------------------------------------------------------
  for iDimension = 1 : nDimensions
    
    % If it is not possible to use a predefined incrementation scheme,
    % incrementation tables are created for each individual dimension, that
    % contain all the changes to event lengths for all events in the
    % corresponding dimension. This also allows to use non-linear
    % increments
    if ~IncrementationScheme
      IncrementationTable = zeros(nEvents,Vary.Points(iDimension));
    end
    
    if nDimensions == 1
      if isfield(Exp,'Dim')
        Field2Get = 'Dim';
      else
        Field2Get = 'Dim1';
      end
    else
      Field2Get = ['Dim' num2str(iDimension)];
    end
    
    % Scans all the lines of the Dim field. Each line can contain multiple
    % events/fields 
    for iLine = 1 : size(Exp.(Field2Get),1)
      % Gets the string, that lists the events/fields that are to be
      % changed
      FullString = Exp.(Field2Get){iLine,1};
      
      % changed values are seperated with commas
      SplitStrings = regexp(FullString,',','split');
      
      % Loops over all entries in the current line
      for iModifiedEvent = 1 : length(SplitStrings)
              
        % modified pulses are defined through p#.Field and if a '.' is
        % found, the string is split again. Delays are ignored this way
        Strings = regexp(SplitStrings{iModifiedEvent},'\.','split');
        
        EventType = Strings{1}(1);
        EventSpecificIndex = str2double(Strings{1}(2:end));
        
        % -----------------------------------------------------------------
        % Different Processing for a pulse 'p' and a free evolution
        % event/delay 'd'
        % -----------------------------------------------------------------
        switch EventType
          case 'p'         
            % Convert the index as provided in the Dimension structure to
            % an eventnumber and pulsenumber
            EventNumber = Pulses{EventSpecificIndex}.EventIndex;
            PulseNumber = EventSpecificIndex;
            
            % Gets the field that is to be modified
            Field = Strings{2};
            
            % Catch if user defines pulse length as 't' instead of 'tp'
            if strcmp(Field,'t')
              Field = 'tp';
            end
            
            switch Field
              % If the field is 'Position', the surrounding events (which
              % have to be delays) are changed in length. This is written
              % to the incrementation table
              case 'Position'
                SurroundingEvents = [EventNumber-1 EventNumber+1];
                
                % get the Increment
                dt = Exp.(Field2Get){iLine,2};
                          
                if ~IncrementationScheme
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % At this point, only a linear increment is processed, in
                  % the future, the ability to use nonlinear increments has
                  % to be added here, by providing a vector with values
                  % instead of a scalar in the Dimension field. This should
                  % only require a few lines of code
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  IncrementationTable(SurroundingEvents(1),:) = IncrementationTable(SurroundingEvents(1),:) + (0:Vary.Points(iDimension)-1)*dt;
                  IncrementationTable(SurroundingEvents(2),:) = IncrementationTable(SurroundingEvents(2),:) - (0:Vary.Points(iDimension)-1)*dt;
                else
                  % write the linear matrix for use with incrementation scheme here
                end
                

              otherwise
                % If not the position is changed it is a pulse parameter.
                % All pulse parameters are first stored in  a seperate
                % structure, called PulseModifications. Each dimension can
                % add fields and values to it. Only after all dimensions
                % have been checked for pulse modifications, the pulses can
                % be calculated and stored in the Vary structure
                if isempty(PulseModifications{PulseNumber})
                  PulseModifications{PulseNumber} = {iDimension Field Exp.(Field2Get){iLine,2}};
                else
                  n = size(PulseModifications{PulseNumber},1);
                  PulseModifications{PulseNumber}{n+1,1} = iDimension;
                  PulseModifications{PulseNumber}{n+1,2} = Field;
                  PulseModifications{PulseNumber}{n+1,3} = Exp.(Field2Get){iLine,2};
                end
            end
            
          case 'd'
            % If a delay is changed, the incrementation/decrementation is
            % written to the incrementation table of the corresponding
            % dimension
            EventNumber = DelayIndices(EventSpecificIndex);
            
            % Increment
            dt = Exp.(Field2Get){iLine,2};
            
            if ~IncrementationScheme
              IncrementationTable(EventNumber(1),:) = IncrementationTable(EventNumber(1),:) + (0:Vary.Points(iDimension)-1)*dt;
            else
              %write the linear matrix here
            end
        end
      end
    end
    
    % Stores the IncrementationTable dimension specific
    if ~IncrementationScheme && any(any(IncrementationTable))
      Vary.IncrementationTable{iDimension} = IncrementationTable;
    else
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Store Directly as Vary.IncrementationTable for Incrementation
      % Schemes
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
  end
  
  % -----------------------------------------------------------------------
  % The following part checks for pulse overlap and precomputes the wave
  % forms for pulses that are variied
  % -----------------------------------------------------------------------
  nDataPoints = prod(Vary.Points);
  DimensionIndices = ones(1,nDimensions);
  
  % For each data point, the originial/starting values for the pulses are
  % required
  InitialPulses = Pulses;
  
  % Each pulse that is modified will have all its waveforms stored in
  % Vary.Pulses{iPulse}. Vary.Pulses{iPulse} is a cell array with the
  % dimensionality of points in the dimension that it is being changed in.
  % If for example the pulse is changed along the first dimension the size
  % of Vary.Pulses{iPulse} is [nDim1 1 1 ...], if it is changed along the
  % second dimension [1 nDim2 1 1 ...] or if along all dimension [nDim1
  % nDim2 nDim3 nDim4 ...]
  Vary.Pulses = cell(1,nPulses);
  
  % Loop over all DataPoints/Aquisitions and check for pulse overlap and
  % compute pulse shapes if they are changed
  for iDataPoint = 1 : nDataPoints
    % Load starting values for pulses and event lengths
    Pulses = InitialPulses;
    EventLengths = Exp.t;
    
    % ---------------------------------------------------------------------
    % First we need to loop over all pulses and check them for 
    % modifications
    % ---------------------------------------------------------------------
    for iPulse = 1 : nPulses
      % if this pulse is changed, set the values for the pulse parameters
      % accordingly
      if ~isempty(PulseModifications{iPulse})
        % Create and index for storing the pulse shapes in
        % Vary.Pules{iPulse}
        Pulses{iPulse}.ArrayIndex = ones(1,length(DimensionIndices));
        for iModification = 1 : size(PulseModifications{iPulse},1)
          % Load Modifications
          Dimension = PulseModifications{iPulse}{iModification,1};
          Field = PulseModifications{iPulse}{iModification,2};
          Increment = PulseModifications{iPulse}{iModification,3};
          % Write modifications to pulse structure
          Pulses{iPulse}.(Field) = Pulses{iPulse}.(Field) + Increment*(DimensionIndices(Dimension)-1);
          % Adapt indexing according to dimension
          Pulses{iPulse}.ArrayIndex(Dimension) = DimensionIndices(Dimension);
        end
        
        % Convert array into cell for indexing
        ArrayIndex = num2cell(Pulses{iPulse}.ArrayIndex);
        
        % Compute Wave form and store it
        for iPCstep = 1 : size(Pulses{iPulse}.PhaseCycle,1)
          Pulses{iPulse}.Phase = Pulses{iPulse}.Phase+Pulses{iPulse}.PhaseCycle(iPCstep,1);
          [t,IQ] = pulse(Pulses{iPulse});
          Vary.Pulses{iPulse}.IQs{ArrayIndex{:}}(iPCstep,:) = IQ;
        end
        Vary.Pulses{iPulse}.ts{ArrayIndex{:}} = t;

        
        % Write pulse length to EventLenghts        
        EventLengths(Pulses{iPulse}.EventIndex) = t(end);
      end
    end
    
    % ---------------------------------------------------------------------
    % Now we loop over all dimensions and load the delays from the
    % IncrementationTables
    %----------------------------------------------------------------------
    for iDimension = 1 : nDimensions
      % Check if IncrementationTable is not empty
      if ~IncrementationScheme && ~isempty(Vary.IncrementationTable{iDimension})
        % Find Events that are modified...
        ModifiedEvents = find(Vary.IncrementationTable{iDimension}(:,DimensionIndices(iDimension)));
        % ... and change them in EventLenghts
        if ~isempty(ModifiedEvents)
          for i = 1 : length(ModifiedEvents)
            EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable{iDimension}(ModifiedEvents(i),DimensionIndices(iDimension));
          end
        end
      end
    end
     
    % Reorder Sequence and check for pulse overlap
    NewSequence = reorder_events(EventLengths,isPulse);
    
    % Assert that if events are being moved in the sequence, the values for
    % Detection and Relaxation are the same
    for iEvent = 1 : nEvents
      if Events{NewSequence(iEvent)}.Detection ~= Events{(iEvent)}.Detection
        MessagePart1 = ['Due to a moving pulse, the events ' num2str(iEvent) ' and ' num2str(NewSequence(iEvent))];
        MessagePart2 = ' are being interchanged, but they do not have the same setting with respect to detection.';
        Message = [MessagePart1 MessagePart2];
        error(Message);
      end
      if Events{NewSequence(iEvent)}.Relaxation ~= Events{(iEvent)}.Relaxation
        MessagePart1 = ['Due to a moving pulse, the events ' num2str(iEvent) ' and ' num2str(NewSequence(iEvent))];
        MessagePart2 = ' are being interchanged, but they do not have the same setting with respect to relaxation.';
        Message = [MessagePart1 MessagePart2];
        error(Message);
      end
    end
    
    % Increment Dimension index
    for d = nDimensions:-1:1
      if DimensionIndices(d) < Vary.Points(d)
        DimensionIndices(d) = DimensionIndices(d)+1;
        break;
      else
        DimensionIndices(d) = 1;
      end
    end
    
  end 
end