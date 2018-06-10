function [Events, Vary, Opt] = sequencer(Exp,Opt)
% This function creates the Event and Vary Structures

% -------------------------------------------------------------------------
% Pre-Processing
% -------------------------------------------------------------------------
Vary = [];

% Check if resonator is 
if isfield(Exp,'Resonator')
  IncludeResonator = true;
  
  if ~isfield(Exp,'mwFreq')
    error('For using a resonator, the field Exp.mwFreq needs to be provided, and Exp.Frequency needs to be defined in relation to that.')
  end
  
  if ~isfield(Exp.Resonator,'nu') && ~isfield(Exp.Resonator,'TransferFunction') && ~isfield(Exp.Resonator,'nu0') && ~isfield(Exp.Resonator,'QL')
    error('In order to use a resonator either nu0 and QL or nu and TransferFunction need to be defined.')
    % Looks for frequency axis nu and transfer function
  elseif (isfield(Exp.Resonator,'nu') && ~isfield(Exp.Resonator,'TransferFunction')) || (~isfield(Exp.Resonator,'nu') && isfield(Exp.Resonator,'TransferFunction'))
    error('Either Exp.Resonator.nu or Exp.Resonator.TransferFunction is missing')
  elseif isfield(Exp.Resonator,'nu') && isfield(Exp.Resonator,'TransferFunction')
    Resonator.Arg1 = Exp.Resonator.nu;
    Resonator.Arg2 = Exp.Resonator.TransferFunction;
  elseif (isfield(Exp.Resonator,'nu0') && ~isfield(Exp.Resonator,'QL')) || (~isfield(Exp.Resonator,'nu0') && isfield(Exp.Resonator,'QL'))
    % Looks for center frequency nu0 and loaded Qualityfactor
    error('Either Exp.Resonator.nu0 or Exp.Resonator.QL is missing')
  elseif isfield(Exp.Resonator,'nu0') && isfield(Exp.Resonator,'QL')
    Resonator.Arg1 = Exp.Resonator.nu0;
    Resonator.Arg2 = Exp.Resonator.QL;
  end
  
  % if no mode for the resonator incorporation is given, 'simulate' is
  % assumed by default
  if isfield(Exp.Resonator,'Mode')
    if any(strcmp(Exp.Resonator.Mode,{'simulate' 'compensate'}))
      Resonator.Arg3 = Exp.Resonator.Mode;
    else
      error('Resonator.Mode must be ''simulate'' or ''compensate''.')
    end
  else
    Resonator.Arg3 = 'simulate';
  end
else
  IncludeResonator = false;
end

FreqShift = 0;

% Set FrequencyShift and store them in Opt
if isfield(Exp,'mwFreq')
  FreqShift = FreqShift + Exp.mwFreq;
end

if isfield(Opt,'SimulationMode') && any(strcmp(Opt.SimulationMode,{'LabFrame' 'FrameShift'}))
  switch Opt.SimulationMode
    case 'LabFrame'
      FreqShift = FreqShift; % Do nothing
      FrameShift = false;
    case 'FrameShift'
      if ~isfield(Opt,'FrameShift')
        error('If you request ''FrameShift'' with Opt.SimulationMode, you need to provide frequency for shifting with Opt.FrameShift.')
      end
      FreqShift = FreqShift - Opt.FrameShift;
      FrameShift = Opt.FrameShift;
  end
else
  % If no simulation mode is specified, 
  MinFreq = min(abs(Exp.Frequency(:)+FreqShift));
  FrameShift = floor(MinFreq-2);
  if FrameShift > 0
    FreqShift = FreqShift - FrameShift;
  end
  if isfield(Opt,'FrameShift')
    warning('You provided the field Opt.FrameShift, but did not set Opt.SimulationMode accordingly. Opt.FrameShift ignored, FrameShift is computed automatically.')
  end	
end

Opt.FrameShift = FrameShift;
Opt.FreqShift = FreqShift;

% Check if TimeStep exists and if it is sufficient or, if none provided,
% compute a new one
MaxFreq = max(abs(Exp.Frequency(:)+FreqShift));
Nyquist = 2*MaxFreq;
MaxTimeStep = 1/Nyquist/1000; %%%%%%% Time Step is in microseconds and Frequencies in GHz

if isfield(Exp,'TimeStep')
  if Exp.TimeStep > MaxTimeStep
    errMsg = ['Your Time Step (Exp.TimeStep) does not fullfill the Nyquist criterium for the pulses you provided. Adapt it to ' num2str(MaxTimeStep) ' or less.'];
    error(errMsg);
  end
else
 Exp.TimeStep = round(MaxTimeStep/2,2,'significant');
 message = ['Exp.TimeStep not provided, assuming ' num2str(Exp.TimeStep) 'us as timestep'];
 warning(message)
end

Exp.Frequency = Exp.Frequency*1000; % GHz -> MHz

% Create an empty cell array for all the events
Events = cell(1,length(Exp.t));
nEvents = length(Exp.t);

% Variables for bookkeeping of pulses and free evolution events
iPulse = 0;
iDelay = 0;

% A vector to quickly identify pulses, required for the reordering if
% pulses cross during the sequence
isPulse = zeros(1,length(Events));
PulseIndices = [];


% Doing some preliminary checks if all parameters have the correct lengths
if length(Exp.Pulses) ~= length(Exp.t)
  error('The lengths of Exp.t and Exp.Pulses do not match.')
end

if isfield(Opt,'Relaxation')
  if length(Opt.Relaxation) ~= 1 && length(Opt.Relaxation) ~= length(Exp.t)
    error('The lengths of Exp.t and Opt.Relaxation do not match. Length of Opt.Relaxation has to be 1 or the same as Exp.t.')
  end
end

if isfield(Exp,'DetEvents')
  if length(Exp.DetEvents) ~= 1 && length(Exp.DetEvents) ~= length(Exp.t)
    error('The lengths of Exp.t and Exp.DetEvents do not match. Length of Exp.DetEvents has to be 1 or the same as Exp.t.')
  end
end

if isfield(Opt,'StateTrajectories')
  if length(Opt.StateTrajectories) ~= 1 && length(Opt.StateTrajectories) ~= length(Exp.t)
    error('The lengths of Exp.t and Opt.StateTrajectories do not match. Length of Opt.StateTrajectories has to be 1 or the same as Exp.t.')
  end
end

% Setting up data structures for the pulses and events
for iEvent = 1 : nEvents
  if isstruct(Exp.Pulses{iEvent})
    iPulse = iPulse + 1;
    PulseIndices(iPulse) = iEvent; %#ok<AGROW>
    isPulse(iEvent) = true;
  else
    iDelay = iDelay + 1;
    DelayIndices(iDelay) = iEvent; %#ok<AGROW>
    isPulse(iEvent) = false; 
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
  if isPulse(iEvent)
    Pulse = [];
    
    % Gets the PhaseCycle for the current Pulse, if none is provided, phase
    % cycling is switched off for this event
    if isfield(Exp,'PhaseCycle') &&  iPulse <= length(Exp.PhaseCycle) && ~isempty(Exp.PhaseCycle{iPulse})
      ThisPhaseCycle = Exp.PhaseCycle{iPulse};
      nPhaseSteps = size(ThisPhaseCycle,1);
    else
      ThisPhaseCycle = 0;
      nPhaseSteps = 1;
    end
    
    if isfield(Exp.Pulses{iEvent},'IQ')
      
      
      if iscell(Exp.Pulses{iEvent}.IQ)
        UserIQ = Exp.Pulses{iEvent}.IQ{1};
      elseif ismatrix(Exp.Pulses{iEvent}.IQ)
        UserIQ = Exp.Pulses{iEvent}.IQ;
      else
        errMsg = ['The data structure provided for event ' num2str(iEvent) ' is not recognized.'];
        error(errMsg);
      end
      
      [d1, d2] = size(UserIQ);
      if d1 ~= nPhaseSteps && d2 ~= nPhaseSteps
        errMsg = ['The dimensionality of the IQ signal provided for event ' num2str(iEvent) ' is not in agreement with the phasecycle for this pulse. For user-defined waveforms the array must contain all IQs.'];
        error(errMsg);
      elseif d2 == nPhaseSteps
        UserIQ =  UserIQ';
      end
      
      if ~isfield(Exp.Pulses{iEvent},'t')
        errMsg = ['An userdefined IQ was used for event ' num2str(iEvent) ', but the time axis is  missing. Please provide either through Pulse.t.' ];
        error(errMsg);
      else
        if iscell(Exp.Pulses{iEvent}.t)
          Pulse.TimeStep = Exp.Pulses{iEvent}.t{1}(2) - Exp.Pulses{iEvent}.t{1}(1);
          tIQ = Exp.Pulses{iEvent}.t{1};
        else
          Pulse.TimeStep = Exp.Pulses{iEvent}.t(2) - Exp.Pulses{iEvent}.t(1);
          tIQ = Exp.Pulses{iEvent}.t;
        end
        Pulse.userIQ.t = Exp.Pulses{iEvent}.t;
      end
                        
      % Shifts IQ of the pulse if necessary...
      for iPhaseStep = 1 : nPhaseSteps
          Opt.dt = Exp.TimeStep;
          if IncludeResonator
            % if resonator is requested, pulses are elongated due to ringing.
            % the duration of ringing is stored in an additional field
            tOrig = t(end);
            [tIQ,currentIQ] = resonator(t,UserIQ(iPhaseStep,:),Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
            Events{iEvent}.Ringing = t(end) - tOrig;
          else
            currentIQ = UserIQ(iPhaseStep,:);
          end
          [t, ShiftedUserIQ(iPhaseStep,:)] = rfmixer(tIQ,currentIQ,-Opt.FrameShift,'IQshift',Opt);
      end
      Events{iEvent}.IQ =  ShiftedUserIQ;
      
      Exp.t(iEvent) = t(end); % updates the Exp.t structure to make sure pulse overlaps are caught (resonator, moving pulses) in case the originally provided value in Exp.t was not correct
      
      Pulse.userIQ.IQ = Exp.Pulses{iEvent}.IQ;
      Pulse.PhaseCycle = ThisPhaseCycle;
    else
      % ---------------------------------------------------------------------
      % Pulse Specific Fields
      % ---------------------------------------------------------------------
      Pulse = Exp.Pulses{iEvent};
      
      Pulse.tp = Exp.t(iEvent);
      
      Pulse.PhaseCycle = ThisPhaseCycle;
      
      % Gets the frequency band from Exp.Frequency. If only one frequency
      % band is provided it is used for all pulses.
      if size(Exp.Frequency,1) == 1
        Pulse.Frequency = Exp.Frequency;
      elseif size(Exp.Frequency,1) < iPulse
        error('The Frequency Band for Pulse No. %d is missing.',iPulse)
      else
        Pulse.Frequency = Exp.Frequency(iPulse,:);
      end
      
      %     if length(Pulse.Frequency) == 2 && Pulse.Frequency(1) == Pulse.Frequency(2)
      %        Pulse.Frequency = Pulse.Frequency(1);
      %     end
      
      % Gets the flip angle
      if isfield(Exp,'Flip') && iPulse <= length(Exp.Flip)
        Pulse.Flip = Exp.Flip(iPulse);
      elseif ~isfield(Pulse,'Qcrit') && ~isfield(Pulse,'Amplitude')
        error('No Flipangle for Pulse No. %d provided.',iPulse)
      end
      
      % Gets the phase for the pulse, if none is provided, the phase is
      % assumed to be 0
      if isfield(Exp,'Phase') && length(Exp.Phase) >= iPulse
        if isfield(Pulse,'Phase')
          Pulse.Phase = Pulse.Phase + Exp.Phase(iPulse);
        else
          Pulse.Phase = Exp.Phase(iPulse);
        end
      elseif ~isfield(Pulse,'Phase')
        Pulse.Phase = 0;
      end
      

      % Get the time step size if available
      Pulse.TimeStep = Exp.TimeStep;
      
      % Loop over the function that creates the PulseShape, as many times at
      % are necessary to calculate all wave forms for the phase cycling
      for iPCstep = 1 : nPhaseSteps
        Pulse.Phase = Pulse.Phase + Pulse.PhaseCycle(iPCstep,1);
        [t,IQ] = pulse(Pulse);
        if IncludeResonator
          % if resonator is requested, pulses are elongated due to ringing.
          % the duration of ringing is stored in an additional field
          tOrig = t(end);
          [t,IQ] = resonator(t,IQ,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
          Events{iEvent}.Ringing = t(end) - tOrig;
        end
        % Shifts IQ of the pulse if necessary...
        if FreqShift ~= 0
          Opt.dt = Exp.TimeStep;
          [t, IQ] = rfmixer(t,IQ,FreqShift,'IQshift',Opt);
        end
        % ... and stores it in the event structure
        Events{iEvent}.IQ(iPCstep,:) = IQ;
      end
      
    end
    
    % Specify Type in Event structure
    Events{iEvent}.type = 'pulse';
    
    % Store the time axis of the pulse in the Event structure    
    Events{iEvent}.t = t;
    
    % Store the PhaseCycle in the Event structure
    Events{iEvent}.PhaseCycle = ThisPhaseCycle;
       
    % Checks if ComplexExcitation is requested for this Pulse, if not
    % specified Complex Excitation is switched off by default - the
    % excitation operator is being built outside of sequencer
    if ~isfield(Opt,'ComplexExcitation') || isempty(Opt.ComplexExcitation)
      Events{iEvent}.ComplexExcitation = false;
    elseif length(Opt.ComplexExcitation) == 1
      Events{iEvent}.ComplexExcitation = Opt.ComplexExcitation;
    elseif iPulse <= length(Opt.ComplexExcitation)
      Events{iEvent}.ComplexExcitation = Opt.ComplexExcitation(iPulse);
    else
      Events{iEvent}.ComplexExcitation = false;
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
    else
      Events{iEvent}.Relaxation = Opt.Relaxation(iEvent);
    end
  end
  
  % Check if detection is provided, if no detection is requested, detection
  % is switched off
  if ~isfield(Exp,'DetEvents') || isempty(Exp.DetEvents)
      Events{iEvent}.Detection = false;
  else
    if length(Exp.DetEvents) == 1
      Events{iEvent}.Detection = Exp.DetEvents;
    else 
      Events{iEvent}.Detection = Exp.DetEvents(iEvent);
    end   
  end
  
  % Check if Density Matrices are to be stored, if not specified, Density
  % Matrices are not stored
  if ~isfield(Opt,'StateTrajectories')
    Events{iEvent}.StateTrajectories = false;
  else
    if length(Opt.StateTrajectories) == 1
      Events{iEvent}.StateTrajectories = Opt.StateTrajectories;
    else
      Events{iEvent}.StateTrajectories = Opt.StateTrajectories(iEvent);
    end
  end
  
  % Store an empty propagation structure, will be overwritten by thyme
  Events{iEvent}.Propagation = [];
  
  % Store the time step, which will be needed in thyme to calculate time
  % axis and propagators
  Events{iEvent}.TimeStep = Exp.TimeStep;
  
  % Keep track of the frameshift
  Events{iEvent}.FrameShift = FrameShift;
    
end

% -------------------------------------------------------------------------
% Checks for overlap of pulses that are subject to ringing
% -------------------------------------------------------------------------
if IncludeResonator
  for iEvent = PulseIndices
    FollowingEvent = iEvent + 1;
    if FollowingEvent <= length(Exp.t) && strcmp(Events{FollowingEvent}.type,'pulse')
      error('When using a resonator, pulses need to be separated by inter pulse delays to accomodate for ringing from the resonator.')
    elseif FollowingEvent <= length(Exp.t)
      ShortenedDelay = Events{FollowingEvent}.t - Events{iEvent}.Ringing;
      if ShortenedDelay < 0
        Msg = ['Event ' num2str(FollowingEvent) ' (a delay) is too short to accomodate ringing of the preceding pulse.'];
        error(Msg);
      end
    end    
  end  
end

% -------------------------------------------------------------------------
% Creates the Vary structure
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switching incrementation scheme to off for now, this will be an option
% later on, and will decide on how the incremenation tables are stored. For
% the incrementationscheme only linear increments can be used, and the data
% structure can therefore be reduced
IncrementationScheme = false;
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
        error('You provided Exp.Dim, please always provide Dim with a number, e.g Exp.Dim1.')
      end
    end
    
    Field2Get = ['Dim' num2str(iDimension)];
    
    if ~isfield(Exp,Field2Get)
      msg = ['You requested a ' num2str(nDimensions) '-dimensional experiment, but Exp.Dim' num2str(iDimension) ' is missing.'];
      error(msg);
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
        FieldIndex = [];      
        % modified pulses are defined through p#.Field and if a '.' is
        % found, the string is split again. Delays are ignored this way
        Strings = regexp(SplitStrings{iModifiedEvent},'\.','split');
        
        EventType = Strings{1}(1);
        EventSpecificIndex = str2double(Strings{1}(2:end));
        
        if length(Strings) == 2
          pars = regexp(Strings{2},'\(|\)','split');
          Strings{2} = pars{1};
          if length(pars) > 1          
            FieldIndex = str2double(pars{2});
          end
        end
        
        if (length(Strings) ~= 2 || ~strcmp(Strings{2},'IQ'))
          if length(Strings) ~= 2 || ~strcmp(Strings{2},'Frequency') || length(Exp.(Field2Get){iLine,2}) == 1 || ~any(size(Exp.(Field2Get){iLine,2},1) == [1 Exp.nPoints(iDimension)-1])
            if length(Exp.(Field2Get){iLine,2}) ~= 1 && length(Exp.(Field2Get){iLine,2}) ~= Exp.nPoints(iDimension)
              message = ['The number of points provided for Dimension ' num2str(iDimension) ' does not match the length of the vector in the Exp.Dim structure.'];
              error(message);
            end
          end
        end
        
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
              if length(Strings) == 1
                message = ['You requested a pulse to be changed in Exp.' (Field2Get) ' but did not specify the field.'];
                error(message)
              end
              
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
                  
                  if any(SurroundingEvents>nEvents) || any(SurroundingEvents>nEvents)
                    error('Moving pulses can not be the first or last event in your Exp structure.')
                  end
                  
                  % get the Increment
                  dt = Exp.(Field2Get){iLine,2};
                  
                  if ~IncrementationScheme
                    if length(dt) == 1
                      IncrementationTable(SurroundingEvents(1),:) = IncrementationTable(SurroundingEvents(1),:) + (0:Vary.Points(iDimension)-1)*dt;
                      IncrementationTable(SurroundingEvents(2),:) = IncrementationTable(SurroundingEvents(2),:) - (0:Vary.Points(iDimension)-1)*dt;
                    else
                      IncrementationTable(SurroundingEvents(1),1:end) = IncrementationTable(SurroundingEvents(1),1:end) + dt;
                      IncrementationTable(SurroundingEvents(2),1:end) = IncrementationTable(SurroundingEvents(2),1:end) - dt;
                    end
                  else
                    Vary.IncrementationTable(SurroundingEvents(1),iDimension) = Vary.IncrementationTable(SurroundingEvents(1),iDimension) + dt;
                    Vary.IncrementationTable(SurroundingEvents(2),iDimension) = Vary.IncrementationTable(SurroundingEvents(2),iDimension) - dt;
                  end
                  
                case 'IQ'
                  if isempty(PulseModifications{PulseNumber})
                    PulseModifications{PulseNumber} = {iDimension Field [] []};
                  else
                    n = size(PulseModifications{PulseNumber},1);
                    PulseModifications{PulseNumber}{n+1,1} = iDimension;
                    PulseModifications{PulseNumber}{n+1,2} = Field;
                  end
                  
                otherwise
                  % If not the position is changed it is a pulse parameter.
                  % All pulse parameters are first stored in  a seperate
                  % structure, called PulseModifications. Each dimension can
                  % add fields and values to it. Only after all dimensions
                  % have been checked for pulse modifications, the pulses can
                  % be calculated and stored in the Vary structure
                  if isempty(PulseModifications{PulseNumber})
                    if strcmp(Field,'Frequency') 
                      % This is necessary to make sure that 'Frequency' is 
                      % set correctly. In the defintion GHz are being used 
                      % which need to be turned into MHz (line 99) for pulse().
                      PulseModifications{PulseNumber} = {iDimension Field Exp.(Field2Get){iLine,2}*1000 FieldIndex};
                    else
                      PulseModifications{PulseNumber} = {iDimension Field Exp.(Field2Get){iLine,2} FieldIndex};
                    end
                  else
                    n = size(PulseModifications{PulseNumber},1);
                    PulseModifications{PulseNumber}{n+1,1} = iDimension;
                    PulseModifications{PulseNumber}{n+1,2} = Field;
                    if strcmp(Field,'Frequency') 
                      % This is necessary to make sure that 'Frequency' is 
                      % set correctly. In the defintion GHz are being used 
                      % which need to be turned into MHz (line 99) for pulse().
                      PulseModifications{PulseNumber}{n+1,3} = Exp.(Field2Get){iLine,2}*1000;
                    else
                      PulseModifications{PulseNumber}{n+1,3} = Exp.(Field2Get){iLine,2};
                    end
                    PulseModifications{PulseNumber}{n+1,4} = FieldIndex;
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
              if length(dt) == 1
                IncrementationTable(EventNumber(1),:) = IncrementationTable(EventNumber(1),:) + (0:Vary.Points(iDimension)-1)*dt;
              else
                IncrementationTable(EventNumber(1),1:end) = IncrementationTable(EventNumber(1),1:end) + dt;
              end
            else
              Vary.IncrementationTable(EventNumber(1),iDimension) = Vary.IncrementationTable(EventNumber(1),iDimension) + dt;
            end
        end
      end
    end
    
    % Stores the IncrementationTable dimension specific
    if ~IncrementationScheme && any(any(IncrementationTable))
      Vary.IncrementationTable{iDimension} = IncrementationTable;
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
        
        IQindex = 0;
        for iModification = 1 : size(PulseModifications{iPulse},1)
          if strcmp(PulseModifications{iPulse}{iModification,2},'IQ')
            IQindex = IQindex + 1;
          end
        end
        
        if ~IQindex
        
          for iModification = 1 : size(PulseModifications{iPulse},1)
            % Load Modifications
            Dimension = PulseModifications{iPulse}{iModification,1};
            Field = PulseModifications{iPulse}{iModification,2};
            Increment = PulseModifications{iPulse}{iModification,3};
            FieldIndex = PulseModifications{iPulse}{iModification,4};
            % Write modifications to pulse structure
            if length(Increment) == 1 || (strcmp(Field,'Frequency') && size(Increment,1) == 1)
              if isempty(FieldIndex)
                Pulses{iPulse}.(Field) = Pulses{iPulse}.(Field) + Increment*(DimensionIndices(Dimension)-1);
              else
                if length(Increment) == 1
                  Pulses{iPulse}.(Field)(FieldIndex) = Pulses{iPulse}.(Field)(FieldIndex) + Increment*(DimensionIndices(Dimension)-1);
                else
                  Pulses{iPulse}.(Field)(FieldIndex) = Pulses{iPulse}.(Field)(FieldIndex) + Increment(DimensionIndices(Dimension));
                end
              end
            else
              if strcmp(Field,'Frequency')
                Pulses{iPulse}.(Field) = Pulses{iPulse}.(Field) + Increment(DimensionIndices(Dimension),:);
              else
                if isempty(FieldIndex)
                  Pulses{iPulse}.(Field) = Pulses{iPulse}.(Field) + Increment(DimensionIndices(Dimension));
                else
                  Pulses{iPulse}.(Field)(FieldIndex) = Pulses{iPulse}.(Field)(FieldIndex) + Increment(DimensionIndices(Dimension));
                end
              end
            end
            % Adapt indexing according to dimension
            Pulses{iPulse}.ArrayIndex(Dimension) = DimensionIndices(Dimension);
          end
          
          % Convert array into cell for indexing
          ArrayIndex = num2cell(Pulses{iPulse}.ArrayIndex);
                  
          % Compute Wave form and store it
          for iPCstep = 1 : size(Pulses{iPulse}.PhaseCycle,1)
            Pulses{iPulse}.Phase = Pulses{iPulse}.Phase+Pulses{iPulse}.PhaseCycle(iPCstep,1);
            [t,IQ] = pulse(Pulses{iPulse});
            if IncludeResonator
              % if a resonator is present, the ringing duration of each pulse
              % needs to be stored in the vary structure too
              tOrig = t(end);
              [t,IQ] = resonator(t,IQ,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
              Vary.Pulses{iPulse}.Ringing(ArrayIndex{:}) = t(end) - tOrig;
            end
            if FreqShift ~= 0
              Opt.dt = Exp.TimeStep;
              [~, IQ] = rfmixer(t,IQ,FreqShift,'IQshift',Opt);
            end
            % ... and stores it in the vary structure
            Vary.Pulses{iPulse}.IQs{ArrayIndex{:}}(iPCstep,:) = IQ;
          end
          Vary.Pulses{iPulse}.ts{ArrayIndex{:}} = t;

        else
         
          if IQindex < size(PulseModifications{iPulse},1)
            errMsg = ['It is not possible to combine changing pulse parameters with user defined IQ, please check your Dim input for Pulse ' (num2str(iPulse)) '.'];
            error(errMsg)
          end
 
          if size(PulseModifications{iPulse},1) == 1
            Dimension = PulseModifications{iPulse}{1,1};
            IndexToLoad = DimensionIndices(Dimension);
            IndexToSave = ones(length(DimensionIndices));
            IndexToSave(Dimension) = DimensionIndices(Dimension);
            ArrayIndex = num2cell(IndexToSave);
            
            UserIQ = Pulses{iPulse}.userIQ.IQ{IndexToLoad};
            
            tIQ = Pulses{iPulse}.userIQ.t{IndexToLoad};
          else
            IndexToLoad = ones(length(DimensionIndices));
            
            for iModification = 1 : size(PulseModifications{iPulse},1)
              Dimension = PulseModifications{iPulse}{iModification,1};
              IndexToLoad(Dimension) = DimensionIndices(Dimension);
            end
            ArrayIndex = num2cell(IndexToLoad);
            
            UserIQ = Pulses{iPulse}.userIQ.IQ{ArrayIndex{:}};
            
            tIQ = Pulses{iPulse}.userIQ.t{ArrayIndex{:}};
            
          end
          
          nPhaseSteps = size(Pulses{iPulse}.PhaseCycle,1);
          [d1, d2] = size(UserIQ);
          
          if d1 ~= nPhaseSteps && d2 ~= nPhaseSteps
            errMsg = ['The dimensionality of the IQ signal provided for event ' num2str(iEvent) ' is not in agreement with the phasecycle for this pulse. For user-defined waveforms the array must contain all IQs.'];
            error(errMsg);
          elseif d2 == nPhaseSteps
            UserIQ =  UserIQ';
          end
         
          
          % Shifts IQ of the pulse if necessary...
          ShiftedUserIQ = [];
          for iPhaseStep = 1 : nPhaseSteps
            Opt.dt = Exp.TimeStep;
            if IncludeResonator
              % if resonator is requested, pulses are elongated due to ringing.
              % the duration of ringing is stored in an additional field
              tOrig = tIQ(end);
              [tIQ,currentIQ] = resonator(tIQ,UserIQ(iPhaseStep,:),Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
              Vary.Pulses{iPulse}.Ringing(ArrayIndex{:}) = tIQ(end) - tOrig;
            else
              currentIQ = UserIQ(iPhaseStep,:);
            end
            [t, ShiftedUserIQ(iPhaseStep,:)] = rfmixer(tIQ,currentIQ,-Opt.FrameShift,'IQshift',Opt);
          end
          
          Vary.Pulses{iPulse}.IQs{ArrayIndex{:}} = ShiftedUserIQ;

          Vary.Pulses{iPulse}.ts{ArrayIndex{:}} = t;
          
        end

        % Write pulse length to EventLenghts
        if IncludeResonator
          EventLengths(Pulses{iPulse}.EventIndex) = t(end) - Vary.Pulses{iPulse}.Ringing(ArrayIndex{:});
        else
          EventLengths(Pulses{iPulse}.EventIndex) = t(end);
        end
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
        % ... and change them in EventLengths
        if ~isempty(ModifiedEvents)
          for i = 1 : length(ModifiedEvents)
            EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable{iDimension}(ModifiedEvents(i),DimensionIndices(iDimension));
          end
        end
      elseif IncrementationScheme
        % Find Events that are modified...
          ModifiedEvents = find(Vary.IncrementationTable(:,iDimension));
          if ~isempty(ModifiedEvents)
            for i = 1 : length(ModifiedEvents)
              EventLengths(ModifiedEvents(i)) = EventLengths(ModifiedEvents(i)) + Vary.IncrementationTable(ModifiedEvents(i),iDimension)*(DimensionIndices(iDimension)-1);
            end
          end
      end
    end
     
    % Reorder Sequence and check for pulse overlap
    [NewSequence, NewEventLengths] = s_reorder_events(EventLengths,isPulse);
    
    % Check if ringing from the resonator causes pulses to overlap, after
    % they have been reorderd
    if IncludeResonator
      for iPulse = 1 : nPulses
        % get position of the current pulse in the reordered sequence
        ThisEvent = find(NewSequence == PulseIndices(iPulse));
        % Get original Event number of the following event and if...
        if ThisEvent == length(NewSequence)
          break
        end
        FollowingEvent = NewSequence(ThisEvent+1);
        if strcmp(Events{FollowingEvent}.type,'pulse')
          %...the following event is a pulse, create an error
          error('When using a resonator, pulses need to be separated by inter pulse delays to accomodate for ringing from the resonator.')
        else
          %...else the duration of the ringing is being loaded...
          if ~isempty(Vary.Pulses{iPulse})
            Ringing = Vary.Pulses{iPulse}.Ringing(ArrayIndex{:});
          else
            Ringing = Events{PulseIndices(iPulse)}.Ringing;
          end
          %...and the following delay is shortened by the correspoding
          % length
          ShortenedDelay = NewEventLengths(ThisEvent+1) - Ringing;
          if ShortenedDelay < 0
            % if the delay is to short and now becomes negative, an error
            % is returned
            Msg = ['The delay ' num2str(FollowingEvent) ' is too short to accomodate for ringing of the preceeding pulse.'];
            error(Msg);
          end
        end
      end
    end

    % Assert that if events are being moved in the sequence, the values for
    % Detection and Relaxation are the same for events that are being
    % interchanged
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