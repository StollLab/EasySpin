% This function parses the pulse sequence specified in the Exp structure and
% additionally using settings from the Opt structure to create the Event and
% Vary structures needed for the simulation by s_thyme.
% For Opt.SimulationMode="fast", it returns the old-syntax Exp for saffron.
%
% This function is called by saffron and spidyan.
%
% Input:
%   Exp - structure with the following fields
%     Sequence
%     DetSequence
%     DetWindow
%     Dim1, Dim2, etc
%     mwFreq
%     DetTimeStep
%     nPoints
%     mwPolarization
%     ResonatorFrequency
%     ResonatorQL
%     FrequencyResponse
%     ResonatorMode
%     PhaseCycle
%   Opt - structure with the following fields
%     SimulationMode - "fast" or "thyme"
%     StateTrajectories - bool or array of bool to store whether densitry
%         matrices should be stored during evolution periods
%     SimFreq
%     IntTimeStep
%     dt
%     Relaxation  - bool or array of bool to turn on relaxation during
%        evolution periods
%
% Output
%   Events - cell array of structures with the following fields
%     type
%     Ringing
%     IQ
%     t
%     PhaseCycle
%     ComplexExcitation
%     Relaxation
%     Detection
%     StateTrajectories
%     Propagation
%     FrameShift
%     TimeStep
%   Vary - structure with the following fields
%     Points
%     IncrementationTable
%     Pulses
%   Opt - strucure with the following fields
%     SinglePointDetection
%     FrameShift
%     %FreqShift
%   Exp_oldSyntax
%     ...


function varargout = s_sequencer(Exp,Opt)

logmsg(1,'-parsing pulse sequence-----------------------------');

% Validate the input and select a propagation engine
predefinedExperiment = isfield(Exp,'Sequence') && ischar(Exp.Sequence);

message = ''; % initialize message that explains the choice of simulation algorithm - 'fast' or 'general'

if ~isfield(Opt,'SimulationMode'), Opt.SimulationMode = ''; end
origSimulationMode = Opt.SimulationMode;


% Check for fields that prevent the fast simulation mode and switch to thyme
if isfield(Exp,'DetWindow')
  Opt.SimulationMode = 'thyme';
  if predefinedExperiment
    message = addtomessage(message,'predefined experiments can not be run in combination with Exp.DetWindow');
  else
    message = addtomessage(message,'Exp.DetWindow does not work with the fast method');
  end
end
if isfield(Exp,'ResonatorFrequency') || isfield(Exp,'ResonatorQL') || isfield(Exp,'FrequencyResponse')
  Opt.SimulationMode = 'thyme';
  message = addtomessage(message,'Exp definition contains resonator specifics');
end

% Creates the Exp structure in the old saffron syntax
if Opt.SimulationMode=="fast"
  
  if predefinedExperiment
    Exp_oldSyntax = [];
  else
    Exp_oldSyntax = struct;
    
    Sequence = Exp.Sequence;
    
    % Get information about pulse sequence
    pulses = cellfun(@isstruct,Sequence);
    pulsePositions = find(pulses);
    delayPositions = find(~pulses);
    nPulses = length(pulsePositions);
    nDelays = length(delayPositions);
    
    if any(diff(pulsePositions)==0)
      % adjacent pulses cannot be processed by saffron
      Opt.SimulationMode = 'thyme';
      message = addtomessage(message,'two or more pulses are not separated by an interpulse delay');
    end
    
    Flip = zeros(1,nPulses);
    tp = zeros(1,nPulses);
    Phase = ones(1,nPulses);
    
    % Loop over pulses, verify and determine Flip/tp/Phase
    for iPulse = 1:nPulses
      
      pulse_ = Sequence{pulsePositions(iPulse)};
      
      % Make sure it's a rectangular pulse
      if isfield(pulse_,'Type') && pulse_.Type~="rectangular" && pulse_.Type~="rectangular/none"
        message = addtomessage(message,'The fast algorithm only supports monochromatic rectangular pulses');
        Opt.SimulationMode = 'thyme';
      end
      
      % Extract pulse flip angle
      if isfield(pulse_,'Flip')
        Flip(iPulse) = pulse_.Flip/(pi/2);
      else
        Opt.SimulationMode = 'thyme';
        message = addtomessage(message,'flip angles must be provided');
      end
      
      % Extract pulse phase (multiple of pi/2)
      if isfield(pulse_,'Phase')
        Phase(iPulse) = mod(Phase(iPulse) + pulse_.Phase/(pi/2),4);
      end
      
      % Get pulse length (for non-ideal pulses)
      if isfield(pulse_,'tp') && pulse_.tp~=0
        tp(iPulse) = pulse_.tp;
      end
      
    end
    
    % Extract inter-pulse delays
    t = zeros(1,nDelays);
    iDelay = 1;
    for p = delayPositions
      t(iDelay) = Sequence{p};
      iDelay = iDelay + 1;
    end
    
    Exp_oldSyntax.t = t;
    Exp_oldSyntax.Flip = Flip;
    Exp_oldSyntax.tp = tp;
    Exp_oldSyntax.Phase = Phase;
        
    if isfield(Exp,'nPoints')
      nDimensions = length(Exp.nPoints);
      
      Inc = zeros(1,nDelays);
      
      if nDimensions > 2
        message = addtomessage(message,'more than 2 indirect dimensions were provided');
        Opt.SimulationMode = 'thyme';
      end
      dt = zeros(1,nDimensions);
      
      % Loop over the indirect dimensions and check what is being changed -
      % the old saffron engine can only increment delays
      for iDimension = 1 : nDimensions
        dimX = ['Dim' num2str(iDimension)];
        
        % Loop over the lines of Exp.DimX, eg: Exp.DimX = {'d1,d2' 0.4; 'd2' 0.3}
        for iLine = 1 : size(Exp.(dimX),1)
          % Gets the string, that lists the events/fields that are to be
          % changed
          FullString = Exp.(dimX){iLine,1};
          
          % changed values are seperated with commas, eg: {'d1,d2' 0.4}
          SplitStrings = regexp(FullString,',','split');
          
          % loop over the individual comma separated increments: {'d1,d2'
          % 0.4} --> {'d1' 0.4} and {'d2' 0.4}
          for iModifiedEvent = 1 : length(SplitStrings)
            
            Strings = regexp(SplitStrings{iModifiedEvent},'\.','split');
            
            EventType = Strings{1}(1);
            EventSpecificIndex = str2double(Strings{1}(2:end));
            
            if strcmp(EventType,'p')
                % Pulse parameters can not be changed
                message = addtomessage(message,['a pulse parameter is changed along indirect dimension no. ' num2str(iDimension)]);
                Opt.SimulationMode = 'thyme';
            else
              Inc(EventSpecificIndex) = iDimension;
            end
            
            if dt(iDimension) == 0
              dt(iDimension) = Exp.(dimX){iLine,2};
            elseif dt(iDimension) ~= Exp.(dimX){iLine,2}
              Opt.SimulationMode = 'thyme';
              message = addtomessage(message,['the increment of dimension no. ' num2str(iDimension) ' is not linear']);
            end
          end
          
        end
      end
      
      Exp_oldSyntax.Inc = Inc;
      Exp_oldSyntax.dt = dt;
      
      Exp_oldSyntax.nPoints = Exp.nPoints;
    end
  end
end

if Opt.SimulationMode=="fast"  
  varargout = {[],[],Opt,Exp_oldSyntax};
  return
end

if ~strcmp(origSimulationMode,Opt.SimulationMode)
  message = ['The ''thyme'' simulation mode has to be used.' '\n' 'The reason for this is: ' message '\n'];
  logmsg(1,message);
end

% ------------------------------------------------------------------------------
% Pre-processing for the thyme method
% ------------------------------------------------------------------------------

% Issue warning about ignored fields in Exp and Opt
saffronSpecificFieldsExp = {'ExciteWidth','Filter'};
saffronSpecificFieldsOpt = {'TimeDomain','Expand','ProductRule',...
  'EndorMethod','nOffsets','lwOffset','logplot','Window','ZeroFillFactor'};
if any(isfield(Exp,saffronSpecificFieldsExp))
  msg = '';
  for i = 1:numel(saffronSpecificFieldsExp)
    if isfield(Exp,saffronSpecificFieldsExp{i})
      msg = addtomessage(msg,saffronSpecificFieldsExp{i});
    end
  end
  msg = ['The following fields in the Exp structure are specific to the fast algorithm and will be ignored: ' msg];
  disp(msg);
end
if any(isfield(Opt,saffronSpecificFieldsOpt))
  msg = '';
  for i = 1: length(saffronSpecificFieldsOpt)
    if isfield(Opt,saffronSpecificFieldsOpt{i})
      msg = addtomessage(msg,saffronSpecificFieldsOpt{i});
    end
  end
  msg = ['The following fields in the Opt structure are specific to the fast algorithm and will be ignored: ' msg];
  disp(msg);
end

Vary = [];

if ~isfield(Exp,'mwFreq')
  error('Exp.mwFreq is required for Opt.SimulationMethod=''thyme''.');
end

Opt.SinglePointDetection = false;

if predefinedExperiment
  % set up Exp structure from predefined experiment for thyme
  Exp = s_predefinedexperiments(Exp);
end

pulsePositions = find(cellfun(@isstruct,Exp.Sequence));

% Check if pulses have a finite length
idealPulses = false;

% and whether their definition requires Exp.mwFreq to be given
mwFreqrequired = false;

msgtp = 'Please provide pulse lengths (Par.tp) for the following pulses:';
msgFrequency = 'Please provide Par.Frequency for the following pulses:';
for Pos = pulsePositions
  PulseNumber = find(pulsePositions==Pos);
  if ~isfield(Exp.Sequence{Pos},'tp') && ~isfield(Exp.Sequence{Pos},'IQ')
    idealPulses = true;
    msgtp = [msgtp ' ' num2str(PulseNumber) ','];
  end
  if ~isfield(Exp.Sequence{Pos},'Frequency') && ~isfield(Exp,'mwFreq') && ~isfield(Exp.Sequence{Pos},'IQ')
    mwFreqrequired = true;
    msgFrequency = [msgFrequency ' ' num2str(PulseNumber) ','];
  end
  if isfield(Exp.Sequence{Pos},'IQ') && ~isfield(Opt,'IntTimeStep')
    error('If you use userprovided IQs for pulses, Opt.IntTimeStep must be given.')
  end
end

if msgtp(end)==',', msgtp(end) = ''; end
if msgFrequency(end)==',', msgFrequency(end) = ''; end

if idealPulses && mwFreqrequired
  errmsg = ['Real pulses are required for the thyme method.\n' msgtp '\n' msgFrequency];
elseif idealPulses
  errmsg = ['Real pulses are required for the thyme method.\n' msgtp];
elseif mwFreqrequired
  errmsg = ['Real pulses are required for the thyme method.\n' msgFrequency];
else
  errmsg = [];
end

if ~isempty(errmsg)
  error(sprintf(errmsg));
end

% Set up detection
%-------------------------------------------------------------------------------
logmsg(1,'  setting up detection');

if isfield(Exp,'DetWindow')
  if isfield(Exp,'DetSequence')
    warning('You provided both Exp.DetWindow and Exp.DetSequence. Remove one of them.');
  end

  % Ensure that detection window does not overlap with a pulse
  if any(Exp.DetWindow<0) && isstruct(Exp.Sequence{end})
    error('You provided Exp.DetWindow with a negative value, but the last element in Exp.Sequence is a pulse. Detection during a pulse is not possible. Please adapt your detection or append an appropriate free evolution event.')
  elseif min(Exp.DetWindow) < 0 && (abs(min(Exp.DetWindow)) > Exp.Sequence{end})
    error('Your detection window is extending beyond the last free evolution into a pulse. Please shorten detection window or adapt the length of the free evolution.');
  end
  
  % Update Exp.Sequence based on Exp.DetWindow
  if isstruct(Exp.Sequence{end})
    % if the last entry in Sequence is a pulse, add a free-evolution period
    if Exp.DetWindow(1)>0
      Exp.Sequence{end+1} = Exp.DetWindow(1);
    end
  else
    % if the last Sequence entry is a free-evolution period, adjust its time
    Exp.Sequence{end} = Exp.Sequence{end} + Exp.DetWindow(1);
  end
  Opt.SinglePointDetection = isscalar(Exp.DetWindow) || Exp.DetWindow(1) == Exp.DetWindow(2);
  if Opt.SinglePointDetection
    Exp.Sequence{end+1} = 0;
    logmsg(1,'  single-point detection');
  else
    Exp.Sequence{end+1} = diff(Exp.DetWindow);
    logmsg(1,'  transient detection');
  end
  
  % Update Exp.DetSequence
  Exp.DetSequence = zeros(1,numel(Exp.Sequence));
  Exp.DetSequence(end) = true;
  
elseif isfield(Exp,'DetSequence')
  % setting up detection in case of Exp.DetSequence
  logmsg(1,'  setting up detection:');
  if ischar(Exp.DetSequence) || isstring(Exp.DetSequence)
    % parsing strings
    if Exp.DetSequence=="last"
      Exp.DetSequence = zeros(1,length(Exp.Sequence));
      Exp.DetSequence(end) = true;
      logmsg(1,'  detection is active during the last element in Exp.Sequence');
    elseif Exp.DetSequence=="all"
      Exp.DetSequence = ones(1,length(Exp.Sequence));
      logmsg(1,'  all elements in Exp.Sequence are detected');
    else
      msg = 'The string you provided in Exp.DetSequence was not recognized.';
      error(msg);
    end
  else
    if ~isscalar(Exp.DetSequence) && length(Exp.DetSequence) ~= length(Exp.Sequence)
      error('The lengths of Exp.Sequence and Exp.DetSequence do not match. Length of Exp.DetSequence has to be 1 or the same as Exp.Sequence.')
    end
    logmsg(1,'  detection set according to Exp.DetSequence');
  end
  
  % identifying single point detection
  if sum(Exp.DetSequence) == 1 && ~isstruct(Exp.Sequence{Exp.DetSequence==1}) && Exp.Sequence{Exp.DetSequence==1} == 0
    Opt.SinglePointDetection = true;
    logmsg(1,'  single point detection');
  end
  
else
  % default if no detection is given:
  logmsg(1,'  no detection specified, detection is active during the entire sequence');
  Exp.DetSequence = ones(1,numel(Exp.Sequence));
end

% Check if resonator is provided
includeResonator = false;
if any([isfield(Exp,'ResonatorFrequency') isfield(Exp,'ResonatorQL') isfield(Exp,'FrequencyResponse')])
  logmsg(1,'  resonator present');
  if isfield(Exp,'ResonatorFrequency') && isfield(Exp,'ResonatorQL')
    logmsg(1,['  resonator frequency: ' num2str(Exp.ResonatorFrequency)]);
    logmsg(1,['  resonator QL: ' num2str(Exp.ResonatorQL)]);
    Resonator.Arg1 = Exp.ResonatorFrequency;
    Resonator.Arg2 = Exp.ResonatorQL;
    includeResonator = true;
  elseif isfield(Exp,'ResonatorFrequency')
    error('Exp.ResonatorFrequency provided, but Exp.ResonatorQL is missing')
  elseif isfield(Exp,'ResonatorQL')
    error('Exp.ResonatorQL provided, but Exp.ResonatorFrequency is missing')
  end
  
  if isfield(Exp,'FrequencyResponse')
    if includeResonator
      logmsg(1,'  also found FrequencyResponse of the resonator, ignoring ResonatorFrequency and QL');
    end
    logmsg(1,'  using Exp.FrequencyResponse to simulate resonator.');
    Resonator.Arg1 = Exp.FrequencyResponse(1,:);
    Resonator.Arg2 = Exp.FrequencyResponse(2,:);
    includeResonator = true;
  end
  
  if isfield(Exp,'ResonatorMode')
    if any(strcmp(Exp.ResonatorMode,{'simulate' 'compensate'}))
      logmsg(1,['  resonator mode: ' Exp.ResonatorMode]);
      Resonator.Arg3 = Exp.ResonatorMode;
    else
      error('Resonator.Mode must be ''simulate'' or ''compensate''.')
    end
  else
    Resonator.Arg3 = 'simulate';
    logmsg(1,'  found no Exp.ResonatorMode, incorporating resonator, but not compensating for it');
  end
end


% Determine simulation frame frequency
%-------------------------------------------------------------------------------
logmsg(1,'  determining simulation frame frequency');

if isfield(Exp,'mwFreq')
  freqShift = Exp.mwFreq;
else
  freqShift = 0;  % GHz
end

if isfield(Opt,'SimFreq')
  % User-provided value for simulation frequency shift
  frameShift = Opt.SimFreq;
  freqShift = freqShift - frameShift;
  
else

  % Determine minimum frequency used by pulses
  minFreq = inf;
  for iEvent = pulsePositions
    if isfield(Exp.Sequence{iEvent},'Frequency')
      pulseFreq = Exp.Sequence{iEvent}.Frequency/1e3;  % MHz -> GHz
    else
      pulseFreq = 0;
    end
    minFreq = min(minFreq,freqShift+min(pulseFreq));
  end
  
  % Have at least 2 GHz difference to the lowest frequency (for the frame shift)
  limFreq = 2;  % GHz
  frameShift = floor(minFreq-limFreq);

  % Only shift down, not up - if lab frame frequencies exist that are < 2 GHz
  if frameShift>0
    freqShift = freqShift - frameShift;
  end
end

% Store in Opt output structure
Opt.FrameShift = frameShift;

% Reporting
logmsg(1,'    simulation frame frequency: %d GHz',frameShift);
if frameShift==0
  logmsg(1,'    simulating in the lab frame');
end


% Determine integration time step
%-------------------------------------------------------------------------------
% Check if IntTimeStep exists and if it is sufficient or, if none provided,
% compute a new one
logmsg(1,'  determining minimal required time step');

% Determine maximum frequency used by pulses (relative to sim freq)
maxFreq = -inf;
for iEvent = pulsePositions
  if isfield(Exp.Sequence{iEvent},'Frequency')
    pulseFreq = Exp.Sequence{iEvent}.Frequency/1e3;  % MHz -> GHz
  else
    pulseFreq = 0;
  end
  maxFreq = max(freqShift+max(pulseFreq),maxFreq);
end

NyquistFreq = 2*maxFreq;
NyquistTime = 1/NyquistFreq/1000; % Time Step is in microseconds and Frequencies in GHz

IntStepFactor = 50; % control the size of the integration time step - 1/x of Nyquist
DetStepFactor = IntStepFactor/2; % gets the scaling factor for the detection time step (approx. 1/2 the size of the Nyquist step), in terms of multiples of integration time step
RefIntTimeStep = round(NyquistTime/IntStepFactor,2,'significant'); % The default integration time step
RefDetTimeStep = DetStepFactor*RefIntTimeStep; % The default dection time step

% Validate time step
if isfield(Opt,'IntTimeStep') && ~isfield(Exp,'DetTimeStep')
  
  if Opt.IntTimeStep > NyquistTime
    errMsg = ['Your integration time step (Opt.IntTimeStep) does not fulfill the Nyquist criterion for the pulses you provided. For best results, change it to ' num2str(RefIntTimeStep, '%10.1e') ' µs or less.'];
    error(errMsg);
  elseif Opt.IntTimeStep > NyquistTime/(1/2*IntStepFactor)
    warnMsg = ['Although your integration time step (Opt.IntTimeStep) fulfills the Nyquist criterion for the pulses you provided, it might not be small enough for accurate results. You might want to change it to ' num2str(RefIntTimeStep, '%10.1e') ' µs or less.'];
    warning(warnMsg);
  end
  
  logmsg(1,'  automatically assuming a suitable detection time step');
  Exp.DetTimeStep = RefDetTimeStep;
  
elseif ~isfield(Opt,'IntTimeStep') && isfield(Exp,'DetTimeStep')
  
  logmsg(1,'  computing an integration time step that fits the provided Exp.DetTimeStep');
  factor = ceil(Exp.DetTimeStep/RefIntTimeStep);
  Opt.IntTimeStep = Exp.DetTimeStep/factor; 
  
elseif isfield(Opt,'IntTimeStep') && isfield(Exp,'DetTimeStep')
  
  if mod(Opt.IntTimeStep,Exp.DetTimeStep)
    errMsg = 'The provided integration time step (Opt.IntTimeStep) and detection time step (Exp.DetTimeStep) do not match. Exp.DetTimeStep has to be a multiple of Opt.IntTimeStep.';
    error(errMsg);
  end
  
  if Opt.IntTimeStep > NyquistTime
    errMsg = ['Your integration time step (Opt.IntTimeStep) does not fulfill the Nyquist criterion for the pulses you provided. For best results adapt it to ' num2str(RefIntTimeStep, '%10.1e') ' µs or less.'];
    error(errMsg);
  elseif Opt.IntTimeStep > NyquistTime/(2/3*IntStepFactor)
    warnMsg = ['Although your integration time step (Opt.IntTimeStep) fulfills the Nyquist criterion for the pulses you provided, it might not be small enough for accurate results. You might want to adapt it to ' num2str(RefIntTimeStep, '%10.1e') ' µs or less.'];
    warning(warnMsg);
  end
 
else
  
  logmsg(1,'  automatically assuming a suitable integration time step');
  Opt.IntTimeStep = RefIntTimeStep;
  Exp.DetTimeStep = RefDetTimeStep;
  
end

logmsg(1,'  the integration time step is %0.2e microseconds',Opt.IntTimeStep);
logmsg(1,'  the detection time step is %0.2e microseconds',Exp.DetTimeStep);

% Create an empty cell array for all the events
Events = cell(1,length(Exp.Sequence));
nEvents = length(Exp.Sequence);

Intervals = zeros(1,length(Exp.Sequence));

% making sure that relaxation is defined globally of for all individual
% elements in Exp.Sequence
if isfield(Opt,'Relaxation')
  if ~isscalar(Opt.Relaxation) && length(Opt.Relaxation) ~= length(Exp.Sequence)
    error('The lengths of Exp.Sequence and Opt.Relaxation do not match. Length of Opt.Relaxation has to be 1 or the same as Exp.Sequence.')
  end
end

% making sure that state trajectories are defined globally of for all
% individual elements in Exp.Sequence
if isfield(Opt,'StateTrajectories')
  if ~isscalar(Opt.StateTrajectories) && length(Opt.StateTrajectories) ~= length(Exp.Sequence)
    error('The lengths of Exp.Sequence and Opt.StateTrajectories do not match. Length of Opt.StateTrajectories has to be 1 or the same as Exp.Sequence.')
  end
end

% Setting up data structures for the pulses and events
logmsg(1,'  parsing Exp.Sequence:');

% Variables for bookkeeping of pulses and free evolution events
% A vector to quickly identify pulses, required for the reordering if
% pulses cross during the sequence
isPulse = cellfun(@isstruct,Exp.Sequence);
pulseIndices = find(isPulse);
delayIndices = find(~isPulse);
nPulses = length(pulseIndices);
nDelays = length(delayIndices);

if isfield(Exp,'DetWindow')
  logmsg(1,'  found %d pulse(s), %d free evolution period(s) and a detection window',nPulses,nDelays-1);
else
  logmsg(1,'  found %d pulse(s) and %d free evolution period(s)',nPulses,nDelays);
end

if isfield(Exp,'mwPolarization')
  if ~ischar(Exp.mwPolarization) && ~isstring(Exp.mwPolarization)
    error('Exp.mwPolarization has to be ''linear'' or ''circular''.');
  end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create the Events structure
% -------------------------------------------------------------------------
if nPulses > 0
  logmsg(1,'  computing waveforms and setting up the event structures');
else
  logmsg(1,'  setting up the event structures');
end

pulses = cell(1,nPulses);
iPulse = 1;
for iEvent = 1:length(Exp.Sequence)
  if isPulse(iEvent)
    Events{iEvent}.type = 'pulse';
    
    % Gets the PhaseCycle for the current Pulse, if none is provided, phase
    % cycling is switched off for this event
    if isfield(Exp,'PhaseCycle') &&  iPulse <= length(Exp.PhaseCycle) && ~isempty(Exp.PhaseCycle{iPulse})
      thisPhaseCycle = Exp.PhaseCycle{iPulse};
      nPhaseSteps = size(thisPhaseCycle,1);
    else
      thisPhaseCycle = 0;
      nPhaseSteps = 1;
    end
    
    if ~isfield(Exp.Sequence{iEvent},'IQ') % Pulse.IQ is being used for userdefined IQs
      % ---------------------------------------------------------------------
      % Pulse specific fields
      % ---------------------------------------------------------------------
      
      % if no pulse type is provided, the default is a rectangular pulse
      if ~isfield(Exp.Sequence{iEvent},'Type')
        
        % First check for field frequency and correct it
        if ~isfield(Exp.Sequence{iEvent},'Frequency')
          % if no frequency is defined, frequency is set to 0 (Exp.mwFreq
          % will be added later)
          Exp.Sequence{iEvent}.Frequency = 0;
          Exp.Sequence{iEvent}.Type = 'rectangular';
          if ~isfield(Exp,'mwFreq')
            error('If you do not provide a frequency for the pulse, you need to give Exp.mwFreq.')
          end
        elseif ~isscalar(Exp.Sequence{iEvent}.Frequency) && (Exp.Sequence{iEvent}.Frequency(2) ~= Exp.Sequence{iEvent}.Frequency(1))
          % make sure that, if a frequency was provided, it is not a frequency sweep
          msg = ['Pulse at position ' num2str(iEvent) ' in Exp.Sequence: no Pulse.Type specified, assuming a monochromatice rectangular pulse, but the field Pulse.Frequency looks like a frequency-swept pulse.'];
          error(msg);
        end
      elseif strcmp(Exp.Sequence{iEvent}.Type,'rectangular') && ~isfield(Exp.Sequence{iEvent},'Frequency')
        Exp.Sequence{iEvent}.Frequency = 0;
      end
      
      Pulse = Exp.Sequence{iEvent};
      
      if ~isfield(Pulse,'tp')
        error('Pulse at position %d in Exp.Sequence: no Pulse.tp specified. Please provide a pulse length.',iEvent);
      end
      
      % Supplement frequency band if it is missing
      %if ~isfield(Pulse,'Frequency')
      %  error('Pulse at position %d in Exp.Sequence: frequency band (Pulse.Frequency) is missing.',iEvent)
      %end
      if ~isfield(Pulse,'Frequency')
        Pulse.Frequency = [0 0];
      end
      
      Pulse.PhaseCycle = thisPhaseCycle;
      
      % Gets the flip angle
      if isfield(Exp.Sequence{iEvent},'Flip')
        Pulse.Flip = Exp.Sequence{iEvent}.Flip;
      elseif ~isfield(Pulse,'Qcrit') && ~isfield(Pulse,'Amplitude')
        error('Flip angle for pulse no. %d is missing.',iPulse)
      end
      
      % Gets the phase for the pulse, if none is provided, the phase is
      % assumed to be 0
      if ~isfield(Pulse,'Phase')
        Pulse.Phase = 0;
      end
      
      % Get the time step
      Pulse.TimeStep = Opt.IntTimeStep;
            
      % Loop over the function that creates the PulseShape, as many times at
      % are necessary to calculate all wave forms for the phase cycling
      phase0 = Pulse.Phase;
      for iPCstep = 1 : nPhaseSteps
        Pulse.Phase = phase0 + Pulse.PhaseCycle(iPCstep,1);
        [t,IQ] = pulse(Pulse);
        if includeResonator
          % if resonator is requested, pulses are elongated due to ringing.
          % the duration of ringing is stored in an additional field
          tOrig = t(end);
          [t,IQ] = resonator(t,IQ,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
          Events{iEvent}.Ringing = t(end) - tOrig;
        end
        % Shifts IQ of the pulse if necessary...
        if freqShift~=0
          Opt.dt = Opt.IntTimeStep;
          [t, IQ] = rfmixer(t,IQ,freqShift,'IQshift',Opt);
        end
        % ... and stores it in the event structure
        Events{iEvent}.IQ(iPCstep,:) = IQ;
      end
      
    else
      % in case user provided their own IQ(s)
      
      % determine input format:
      if iscell(Exp.Sequence{iEvent}.IQ)
        userIQ = Exp.Sequence{iEvent}.IQ{1};
      elseif ismatrix(Exp.Sequence{iEvent}.IQ)
        userIQ = Exp.Sequence{iEvent}.IQ;
      else
        errMsg = ['The data structure of the userdefined IQ of pulse on position ' num2str(iEvent) ' in Exp.Sequence is not recognized.'];
        error(errMsg);
      end
      
      [d1, d2] = size(userIQ);
      % look for phase cycle and if found, verify that the IQ contains the
      % phase cycle
      if d1 ~= nPhaseSteps && d2 ~= nPhaseSteps
        errMsg = ['The dimensionality of the IQ signal provided of the pulse on position ' num2str(iEvent) ' in Exp.Sequence  is not in agreement with the phasecycle for this pulse. For user-defined waveforms the array must contain all IQs.'];
        error(errMsg);
      elseif d2 == nPhaseSteps
        userIQ =  userIQ';
      end
      
      % look for the time axis Pulse.t that corresponds to Pulse.IQ
      if ~isfield(Exp.Sequence{iEvent},'t')
        errMsg = ['A userdefined IQ was used for the pulse on position ' num2str(iEvent) ' in Exp.Sequence but the time axis is missing. Please provide it through Pulse.t.' ];
        error(errMsg);
      else
        if iscell(Exp.Sequence{iEvent}.t)
          Pulse.TimeStep = Exp.Sequence{iEvent}.t{1}(2) - Exp.Sequence{iEvent}.t{1}(1);
          tIQ = Exp.Sequence{iEvent}.t{1};
        else
          Pulse.TimeStep = Exp.Sequence{iEvent}.t(2) - Exp.Sequence{iEvent}.t(1);
          tIQ = Exp.Sequence{iEvent}.t;
        end
        Pulse.userIQ.t = Exp.Sequence{iEvent}.t;
      end
      
      % Shifts IQ of the pulse if necessary...
      for iPhaseStep = 1 : nPhaseSteps
        Opt.dt = Opt.IntTimeStep;
        if includeResonator
          % but first if resonator is requested, pulses are elongated due to ringing.
          % the duration of ringing is stored in an additional field
          tOrig = t(end);
          [tIQ,currentIQ] = resonator(t,userIQ(iPhaseStep,:),Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
          Events{iEvent}.Ringing = t(end) - tOrig;
        else
          currentIQ = userIQ(iPhaseStep,:);
        end
        [t, shiftedUserIQ(iPhaseStep,:)] = rfmixer(tIQ,currentIQ,-Opt.FrameShift,'IQshift',Opt);
      end
      Events{iEvent}.IQ =  shiftedUserIQ;
      
      Pulse.userIQ.IQ = Exp.Sequence{iEvent}.IQ;
      Pulse.PhaseCycle = thisPhaseCycle;
      
    end
    
    % Store the time axis of the pulse in the Events structure
    Events{iEvent}.t = t;
    
    % Store the PhaseCycle in the Event structure
    Events{iEvent}.PhaseCycle = thisPhaseCycle;
    
    % Checks if ComplexExcitation is requested for this pulse, if not
    % specified ComplexExcitation is switched off by default - the
    % excitation operator is being built outside of s_sequencer
    if isfield(Exp,'mwPolarization') && ~isempty(Exp.mwPolarization) && strcmp(Exp.mwPolarization,'circular')
      Events{iEvent}.ComplexExcitation = true;
    else
      Events{iEvent}.ComplexExcitation = false;
    end
    
    % Temporarily store pulse paramaters to avoid reassigning them for creating the
    % vary table
    Pulse.EventIndex = iEvent;
    pulses{iPulse} = Pulse;
    iPulse = iPulse + 1;
    
    Intervals(iEvent) = t(end);
  else
    % ---------------------------------------------------------------------
    % Fields specific for delays/free
    % ---------------------------------------------------------------------
    Events{iEvent}.type = 'free evolution';
    Events{iEvent}.t = Exp.Sequence{iEvent};
    
    Intervals(iEvent) = Exp.Sequence{iEvent};
  end
  
  % -----------------------------------------------------------------------
  % General fields
  % -----------------------------------------------------------------------
  % The following fields need to be defined for pulses and free evolution events
  
  % Check if Relaxation is requested for Events, otherwise turn it off
  if ~isfield(Opt,'Relaxation')
    Events{iEvent}.Relaxation = false;
  else
    if isscalar(Opt.Relaxation)
      Events{iEvent}.Relaxation = Opt.Relaxation;
    else
      Events{iEvent}.Relaxation = Opt.Relaxation(iEvent);
    end
  end
  
  % Check if detection is requested, otherwise turn it off
  if ~isfield(Exp,'DetSequence') || isempty(Exp.DetSequence)
    Events{iEvent}.Detection = false;
  else
    if isscalar(Exp.DetSequence)
      Events{iEvent}.Detection = Exp.DetSequence;
    else
      Events{iEvent}.Detection = Exp.DetSequence(iEvent);
    end
  end
  
  % Check if density matrices are requested, otherwise turn off
  if ~isfield(Opt,'StateTrajectories')
    Events{iEvent}.StateTrajectories = false;
  else
    if isscalar(Opt.StateTrajectories)
      Events{iEvent}.StateTrajectories = Opt.StateTrajectories;
    else
      Events{iEvent}.StateTrajectories = Opt.StateTrajectories(iEvent);
    end
  end
  
  % Store an empty propagation structure, will be overwritten by thyme
  Events{iEvent}.Propagation = [];
  
  % Store the time step, which will be needed in thyme to calculate time
  % axis and propagators
  Events{iEvent}.TimeStep = Opt.IntTimeStep;
  
  % Keep track of the frameshift
  Events{iEvent}.FrameShift = frameShift;
  
end

% -------------------------------------------------------------------------
% Checks for overlap of pulses that are subject to ringing
% -------------------------------------------------------------------------
if includeResonator
  logmsg(1,'  checking for pulse overlap due to ringing from resonator:');
  
  for iEvent = pulseIndices
    followingEvent = iEvent + 1;
    if followingEvent <= length(Exp.Sequence) && strcmp(Events{followingEvent}.type,'pulse')
      error('When using a resonator, pulses need to be separated by inter pulse delays to accommodate for ringing from the resonator.')
    elseif followingEvent <= length(Exp.Sequence)
      ShortenedDelay = Events{followingEvent}.t - Events{iEvent}.Ringing;
      if ShortenedDelay < 0
        Msg = ['Event ' num2str(followingEvent) ' (a delay) is too short to accommodate ringing of the preceding pulse.'];
        error(Msg);
      end
    end
  end
  logmsg(1,'  all good!');
end

% -------------------------------------------------------------------------
% Creates the Vary structure
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switching incrementation scheme to off for now, this will be an option
% later on, and will decide on how the incremenation tables are stored. For
% the incrementationscheme only linear increments can be used, and the data
% structure can therefore be reduced.
incrementationScheme = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(Exp,'nPoints')
  logmsg(1,'-validating indirect dimensions------------------------');
  nDimensions = length(Exp.nPoints);
  Vary.Points = Exp.nPoints;
  
  % Layout the Incrementation data structure
  if incrementationScheme
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
  logmsg(1,'  found %d indirect dimension(s)',nDimensions);
  logmsg(1,'  with a total %d data points',prod(Vary.Points));
  for iDimension = 1 : nDimensions
    logmsg(1,'  parsing dimension no. %d',iDimension);
    % If it is not possible to use a predefined incrementation scheme,
    % incrementation tables are created for each individual dimension, that
    % contain all the changes to event lengths for all events in the
    % corresponding dimension. This also allows to use non-linear
    % increments
    if ~incrementationScheme
      IncrementationTable = zeros(nEvents,Vary.Points(iDimension));
    end
    
    if nDimensions == 1
      if isfield(Exp,'Dim')
        error('You provided Exp.Dim, please always provide Dim with a number, e.g Exp.Dim1.')
      end
    end
    
    dimX = ['Dim' num2str(iDimension)];
    
    if ~isfield(Exp,dimX)
      msg = ['You requested a ' num2str(nDimensions) '-dimensional experiment, but Exp.Dim' num2str(iDimension) ' is missing.'];
      error(msg);
    end
    
    % Scans all the lines of the Dim field. Each line can contain multiple
    % events/fields
    for iLine = 1 : size(Exp.(dimX),1)
      % Gets the string, that lists the events/fields that are to be
      % changed
      FullString = Exp.(dimX){iLine,1};
      
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
        
        % The field index allows to increment a value of a pulse parameter
        % that is given as vector, eg the inital and final order of a HS
        % pulse: Exp.Dim1 = {'p2.n(2)' 0.5}; --> Field Index gets the (2)
        if length(Strings) == 2
          pars = regexp(Strings{2},'\(|\)','split');
          Strings{2} = pars{1};
          if length(pars) > 1
            FieldIndex = str2double(pars{2});
          end
        end
        
        % check if user provided incrementation vector has the correct
        % length
        if (length(Strings) ~= 2 || ~strcmp(Strings{2},'IQ'))
          % determined whether user wants to modify Pulse.IQ
          if length(Strings) ~= 2 || ~strcmp(Strings{2},'Frequency') || isscalar(Exp.(dimX){iLine,2}) || ~any(size(Exp.(dimX){iLine,2},1) == [1 Exp.nPoints(iDimension)-1])
            % determined if field was p1.Frequency
            if ~isscalar(Exp.(dimX){iLine,2}) && length(Exp.(dimX){iLine,2}) ~= Exp.nPoints(iDimension)
              % checked for correct length
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
            EventNumber = pulses{EventSpecificIndex}.EventIndex;
            PulseNumber = EventSpecificIndex;
            
            % Gets the field that is to be modified
            if isscalar(Strings)
              message = ['You requested a pulse to be changed in Exp.' (dimX) ' but did not specify the field.'];
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
                surroundingEvents = [EventNumber-1 EventNumber+1];
                
                if any(surroundingEvents>nEvents)
                  error('Moving pulses can not be the first or last event in your Exp structure.')
                end
                
                % get the Increment
                dt = Exp.(dimX){iLine,2};
                
                if ~incrementationScheme
                  if isscalar(dt)
                    IncrementationTable(surroundingEvents(1),:) = IncrementationTable(surroundingEvents(1),:) + (0:Vary.Points(iDimension)-1)*dt;
                    IncrementationTable(surroundingEvents(2),:) = IncrementationTable(surroundingEvents(2),:) - (0:Vary.Points(iDimension)-1)*dt;
                  else
                    IncrementationTable(surroundingEvents(1),1:end) = IncrementationTable(surroundingEvents(1),1:end) + dt;
                    IncrementationTable(surroundingEvents(2),1:end) = IncrementationTable(surroundingEvents(2),1:end) - dt;
                  end
                else
                  Vary.IncrementationTable(surroundingEvents(1),iDimension) = Vary.IncrementationTable(surroundingEvents(1),iDimension) + dt;
                  Vary.IncrementationTable(surroundingEvents(2),iDimension) = Vary.IncrementationTable(surroundingEvents(2),iDimension) - dt;
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
                  PulseModifications{PulseNumber} = {iDimension Field Exp.(dimX){iLine,2} FieldIndex};
                else
                  n = size(PulseModifications{PulseNumber},1);
                  PulseModifications{PulseNumber}{n+1,1} = iDimension;
                  PulseModifications{PulseNumber}{n+1,2} = Field;
                  PulseModifications{PulseNumber}{n+1,3} = Exp.(dimX){iLine,2};
                  PulseModifications{PulseNumber}{n+1,4} = FieldIndex;
                end
            end
            
          case 'd'
            % If a delay is changed, the incrementation/decrementation is
            % written to the incrementation table of the corresponding
            % dimension
            EventNumber = delayIndices(EventSpecificIndex);
            
            % Increment
            dt = Exp.(dimX){iLine,2};
            
            if ~incrementationScheme
              if isscalar(dt)
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
    if ~incrementationScheme && any(any(IncrementationTable))
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
  InitialPulses = pulses;
  
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
  logmsg(1,'  creating the Vary structure that contains all required wave forms and delay changes');
  for iDataPoint = 1 : nDataPoints
    % Load starting values for pulses and event lengths
    pulses = InitialPulses;
    EventLengths = Intervals;
    
    % ---------------------------------------------------------------------
    % First we need to loop over all pulses and check them for
    % modifications
    % ---------------------------------------------------------------------
    for iPulse = 1 : nPulses
      % if this pulse is changed, set the values for the pulse parameters
      % accordingly
      if ~isempty(PulseModifications{iPulse})
        % Create and index for storing the pulse shapes in
        % Vary.Pulses{iPulse}
        pulses{iPulse}.ArrayIndex = ones(1,length(DimensionIndices));
        
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
            if isscalar(Increment) || (strcmp(Field,'Frequency') && size(Increment,1) == 1)
              if isempty(FieldIndex)
                pulses{iPulse}.(Field) = pulses{iPulse}.(Field) + Increment*(DimensionIndices(Dimension)-1);
              else
                if isscalar(Increment)
                  pulses{iPulse}.(Field)(FieldIndex) = pulses{iPulse}.(Field)(FieldIndex) + Increment*(DimensionIndices(Dimension)-1);
                else
                  pulses{iPulse}.(Field)(FieldIndex) = pulses{iPulse}.(Field)(FieldIndex) + Increment(DimensionIndices(Dimension));
                end
              end
            else
              if strcmp(Field,'Frequency')
                pulses{iPulse}.(Field) = pulses{iPulse}.(Field) + Increment(DimensionIndices(Dimension),:);
              else
                if isempty(FieldIndex)
                  pulses{iPulse}.(Field) = pulses{iPulse}.(Field) + Increment(DimensionIndices(Dimension));
                else
                  pulses{iPulse}.(Field)(FieldIndex) = pulses{iPulse}.(Field)(FieldIndex) + Increment(DimensionIndices(Dimension));
                end
              end
            end
            % Adapt indexing according to dimension
            pulses{iPulse}.ArrayIndex(Dimension) = DimensionIndices(Dimension);
          end
          
          % Convert array into cell for indexing
          ArrayIndex = num2cell(pulses{iPulse}.ArrayIndex);
              
          % Compute wave form and store it
          phase0 = pulses{iPulse}.Phase;
          for iPCstep = 1 : size(pulses{iPulse}.PhaseCycle,1)
            pulses{iPulse}.Phase = phase0 + pulses{iPulse}.PhaseCycle(iPCstep,1);
            [t,IQ] = pulse(pulses{iPulse});
            if includeResonator
              % if a resonator is present, the ringing duration of each pulse
              % needs to be stored in the vary structure too
              tOrig = t(end);
              [t,IQ] = resonator(t,IQ,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
              Vary.Pulses{iPulse}.Ringing(ArrayIndex{:}) = t(end) - tOrig;
            end
            if freqShift~=0
              Opt.dt = Opt.IntTimeStep;
              [~, IQ] = rfmixer(t,IQ,freqShift,'IQshift',Opt);
            end
            % ... and stores it in the vary structure
            Vary.Pulses{iPulse}.IQs{ArrayIndex{:}}(iPCstep,:) = IQ;
          end
          
          Vary.Pulses{iPulse}.ts{ArrayIndex{:}} = t;
          
          
        else
          % if user wants to change p1.IQ (user provided IQ)
          
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
            
            userIQ = pulses{iPulse}.userIQ.IQ{IndexToLoad};
            
            tIQ = pulses{iPulse}.userIQ.t{IndexToLoad};
          else
            IndexToLoad = ones(length(DimensionIndices));
            
            for iModification = 1 : size(PulseModifications{iPulse},1)
              Dimension = PulseModifications{iPulse}{iModification,1};
              IndexToLoad(Dimension) = DimensionIndices(Dimension);
            end
            ArrayIndex = num2cell(IndexToLoad);
            
            userIQ = pulses{iPulse}.userIQ.IQ{ArrayIndex{:}};
            
            tIQ = pulses{iPulse}.userIQ.t{ArrayIndex{:}};
            
          end
          
          % get phase cycle
          nPhaseSteps = size(pulses{iPulse}.PhaseCycle,1);
          [d1, d2] = size(userIQ);
          
          % validate phase cycle
          if d1 ~= nPhaseSteps && d2 ~= nPhaseSteps
            errMsg = ['The dimensionality of the IQ signal provided for event ' num2str(iEvent) ' is not in agreement with the phasecycle for this pulse. For user-defined waveforms the array must contain all IQs.'];
            error(errMsg);
          elseif d2 == nPhaseSteps
            userIQ =  userIQ';
          end
          
          % Shifts IQ of the pulse if necessary...
          shiftedUserIQ = [];
          for iPhaseStep = 1 : nPhaseSteps
            Opt.dt = Opt.IntTimeStep;
            if includeResonator
              % if resonator is requested, pulses are elongated due to ringing.
              % the duration of ringing is stored in an additional field
              tOrig = tIQ(end);
              [tIQ,currentIQ] = resonator(tIQ,userIQ(iPhaseStep,:),Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
              Vary.Pulses{iPulse}.Ringing(ArrayIndex{:}) = tIQ(end) - tOrig;
            else
              currentIQ = userIQ(iPhaseStep,:);
            end
            [t, shiftedUserIQ(iPhaseStep,:)] = rfmixer(tIQ,currentIQ,-Opt.FrameShift,'IQshift',Opt);
          end
          
          Vary.Pulses{iPulse}.IQs{ArrayIndex{:}} = shiftedUserIQ;
          
          Vary.Pulses{iPulse}.ts{ArrayIndex{:}} = t;
          
        end
        
        % Write pulse length to EventLengths
        if includeResonator
          EventLengths(pulses{iPulse}.EventIndex) = t(end) - Vary.Pulses{iPulse}.Ringing(ArrayIndex{:});
        else
          EventLengths(pulses{iPulse}.EventIndex) = t(end);
        end
      end
    end
    
    % ---------------------------------------------------------------------
    % Now we loop over all dimensions and load the delays from the
    % IncrementationTables
    %----------------------------------------------------------------------
    
    for iDimension = 1 : nDimensions
      % Check if IncrementationTable is not empty
      if ~incrementationScheme && ~isempty(Vary.IncrementationTable{iDimension})
        % Find Events that are modified...
        modifiedEvents = find(Vary.IncrementationTable{iDimension}(:,DimensionIndices(iDimension)));
        % ... and change them in EventLengths
        if ~isempty(modifiedEvents)
          for i = 1 : length(modifiedEvents)
            EventLengths(modifiedEvents(i)) = EventLengths(modifiedEvents(i)) + Vary.IncrementationTable{iDimension}(modifiedEvents(i),DimensionIndices(iDimension));
          end
        end
      elseif incrementationScheme
        % Find Events that are modified...
        modifiedEvents = find(Vary.IncrementationTable(:,iDimension));
        if ~isempty(modifiedEvents)
          for i = 1 : length(modifiedEvents)
            EventLengths(modifiedEvents(i)) = EventLengths(modifiedEvents(i)) + Vary.IncrementationTable(modifiedEvents(i),iDimension)*(DimensionIndices(iDimension)-1);
          end
        end
      end
    end
    
    % Reorder Sequence and check for pulse overlap
    [newSequence, newEventLengths] = s_reorder_events(EventLengths,isPulse);
    
    % Check if ringing from the resonator causes pulses to overlap, after
    % they have been reorderd
    if includeResonator
      for iPulse = 1:nPulses
        % get position of the current pulse in the reordered sequence
        thisEvent = find(newSequence == pulseIndices(iPulse));
        % Get original Event number of the following event and if...
        if thisEvent == length(newSequence)
          break
        end
        followingEvent = newSequence(thisEvent+1);
        if strcmp(Events{followingEvent}.type,'pulse')
          %...the following event is a pulse, create an error
          error('When using a resonator, pulses need to be separated by inter pulse delays to accommodate for ringing from the resonator.')
        else
          %...else the duration of the ringing is being loaded...
          if ~isempty(Vary.Pulses{iPulse})
            Ringing = Vary.Pulses{iPulse}.Ringing(ArrayIndex{:});
          else
            Ringing = Events{pulseIndices(iPulse)}.Ringing;
          end
          %...and the following delay is shortened by the correspoding
          % length
          ShortenedDelay = newEventLengths(thisEvent+1) - Ringing;
          if ShortenedDelay < 0
            % if the delay is to short and now becomes negative, an error
            % is returned
            Msg = ['The delay ' num2str(followingEvent) ' is too short to accommodate for ringing of the preceeding pulse.'];
            error(Msg);
          end
        end
      end
    end
    
    % Assert that if events are being moved in the sequence, the values for
    % Detection and Relaxation are the same for events that are being interchanged
    for iEvent = 1:nEvents
      if Events{newSequence(iEvent)}.Detection ~= Events{(iEvent)}.Detection
        messagePart1 = ['Due to a moving pulse, the events ' num2str(iEvent) ' and ' num2str(newSequence(iEvent))];
        messagePart2 = ' are being interchanged, but they do not have the same setting with respect to detection.';
        message = [messagePart1 messagePart2];
        error(message);
      end
      if Events{newSequence(iEvent)}.Relaxation ~= Events{(iEvent)}.Relaxation
        messagePart1 = ['Due to a moving pulse, the events ' num2str(iEvent) ' and ' num2str(newSequence(iEvent))];
        messagePart2 = ' are being interchanged, but they do not have the same setting with respect to relaxation.';
        message = [messagePart1 messagePart2];
        error(message);
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

varargout = {Events,Vary,Opt,[]};

logmsg(1,'  pulse sequence parsed successfully!');

end

function newstring = addtomessage(oldstring,toadd)
if isempty(oldstring)
  newstring = toadd;
else
  newstring = [oldstring ', ' toadd];
end
end
