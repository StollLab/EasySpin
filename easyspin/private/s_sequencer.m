function [varargout] = s_sequencer(Exp,Opt)
% This function creates the Event and Vary Structures

logmsg(1,'-validating pulse sequence-----------------------------');

% Validate the input and select a propagation engine
predefinedExperiment = isfield(Exp,'Sequence') && (ischar(Exp.Sequence) || isstring(Exp.Sequence));

message = []; % initialize message that explains the choice of simulation algorithm - 'fast' or 'general'

OrigSimulationMode = Opt.SimulationMode;

saffronSpecificFieldsExp = {'ExciteWidth','Filter'};
saffronSpecificFieldsOpt = {'TimeDomain','Expand','ProductRule',...
  'EndorMethod','nOffsets','lwOffset','logplot','Window','ZeroFillFactor'};
generalFields = {'mwFreq','Field','CrystalOrientation','CrystalSymmetry'};

% Check a few general fields if the fast algorithm can be run at all
if isfield(Exp,'DetWindow')
    Opt.SimulationMode = 'thyme';
    message = addtomessage(message,'Exp.DetWindow does not work with the fast method');
end

if isfield(Exp,'ResonatorFrequency') || isfield(Exp,'ResonatorQL') || isfield(Exp,'FrequencyResponse')
  Opt.SimulationMode = 'thyme';
  message = addtomessage(message,'Exp defintion contains resonator fields');
end

%
if predefinedExperiment && isfield(Exp,'DetWindow')
  Opt.SimulationMode = 'thyme';
  message = addtomessage(message,'predefined experiments can not be run in combination with Exp.DetWindow');
end

% Creates the Exp structure in the old saffron syntax
if ~isfield(Opt,'SimulationMode') || strcmp(Opt.SimulationMode,'fast')
  
  if predefinedExperiment
    Exp_oldSyntax = Exp;
  else
    Exp_oldSyntax = [];
    
    Sequence = Exp.Sequence;
    
    % getting some basic knowledge about the experiment
    Pulses = cellfun(@isstruct,Sequence);
    PulsePositions = find(Pulses);
    DelayPositions = find(~Pulses);
    
    nPulses = length(PulsePositions);
    nDelays = length(DelayPositions);
    
    if any(diff(PulsePositions)==0)
      % neigbouring pulses, can not be processed by saffron
      Opt.SimulationMode = 'thyme';
      message = addtomessage(message,'two or more pulses are not separated by an interpulse delay');
    end
    
    Flip = zeros(1,nPulses);
    tp = zeros(1,nPulses);
    Phase = ones(1,nPulses);
    
    % Loop over pulses, verify and write into saffron specific fields
    for iPulse = 1:nPulses
      
      pulse_ = Sequence{PulsePositions(iPulse)};
      
      % make sure its a rectangular pulse
      if isfield(pulse_,'Type') && ~strcmp(pulse_,'rectangular') && ~strcmp(pulse_,'rectangular/none')
        message = addtomessage(message,'the fast algorithm only supports ideal or monochromatic rectangular pulses');
        Opt.SimulationMode = 'thyme';
      end
      
      % determine flip angle
      if isfield(pulse_,'Flip')
        Flip(iPulse) = pulse_.Flip/(pi/2);
      else
        Opt.SimulationMode = 'thyme';
        message = addtomessage(message,'flip angles must be provided');
      end
      
      if isfield(pulse_,'Phase')
        Phase(iPulse) = mod(Phase(iPulse) + pulse_.Phase/(pi/2),4);
      end
      
      % Get time length for non-ideal pulses
      if isfield(pulse_,'tp') && pulse_.tp ~= 0
        tp(iPulse) = pulse_.tp;
      end
      
    end
    
    % set up saffron fields for delays
    t = zeros(1,nDelays);
    iDelay = 1;
    
    for Pos = DelayPositions
      t(iDelay) = Sequence{Pos};
      iDelay = iDelay + 1;
    end
    
    Exp_oldSyntax.t = t;
    Exp_oldSyntax.Flip = Flip;
    Exp_oldSyntax.tp = tp;
    Exp_oldSyntax.Phase = Phase;
    
    % Populate NewExp with all the other saffron specific fields
    % compares the list of saffron specific fields with
    for iField = 1 : length(generalFields)
      if isfield(Exp,generalFields{iField})
        Exp_oldSyntax.(generalFields{iField}) = Exp.(generalFields{iField});
      end
    end
    
    if isfield(Exp,'nPoints')
      nDimensions = length(Exp.nPoints);
      
      Inc = zeros(1,nDelays);
      
      if nDimensions > 2
        message = addtomessage(message,'more than 2 indirect dimensions were provided');
        Opt.SimulationMode = 'thyme';
      end
      dt = zeros(1,nDimensions);
      
      % loop over the indirect dimensions and check what is being changed -
      % the old saffron engine can only increment delays
      for iDimension = 1 : nDimensions
        Field2Get = ['Dim' num2str(iDimension)];
        
        % Loop over the lines of Exp.DimX, eg: Exp.DimX = {'d1,d2' 0.4; 'd2' 0.3}
        for iLine = 1 : size(Exp.(Field2Get),1)
          % Gets the string, that lists the events/fields that are to be
          % changed
          FullString = Exp.(Field2Get){iLine,1};
          
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
              dt(iDimension) = Exp.(Field2Get){iLine,2};
            elseif dt(iDimension) ~= Exp.(Field2Get){iLine,2}
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

if ~strcmp(OrigSimulationMode,Opt.SimulationMode)
  message = ['The ''thyme'' simulation mode has to be used.' '\n' 'The reason for this was: ' message '\n'];
  logmsg(1,message);
end

if strcmp(Opt.SimulationMode,'thyme')
  if any(isfield(Exp,saffronSpecificFieldsExp))
    msg = [];
    for i = 1: length(saffronSpecificFieldsExp)
      if isfield(Exp,saffronSpecificFieldsExp{i})
        msg = addtomessage(msg,saffronSpecificFieldsExp{i});
      end
    end
    msg = ['The following fields in the Exp structure are specific to the fast algorithm and will be ignored: ' msg];
    disp(msg)
  end 
  if any(isfield(Opt,saffronSpecificFieldsOpt))
    msg = [];
    for i = 1: length(saffronSpecificFieldsOpt)
      if isfield(Opt,saffronSpecificFieldsOpt{i})
        msg = addtomessage(msg,saffronSpecificFieldsOpt{i});
      end
    end
    msg = ['The following fields in the Opt structure are specific to the fast algorithm and will be ignored: ' msg];
    disp(msg)
  end
end

if strcmp(Opt.SimulationMode,'fast')
   
  Exp_oldSyntax.Processed = true;
  varargout{1} = [];
  varargout{2} = [];
  varargout{3} = Opt;
  varargout{4} = Exp_oldSyntax;
  
  return
end

% -------------------------------------------------------------------------
% Pre-Processing for the thyme method
% -------------------------------------------------------------------------

Vary = [];

if ~isfield(Exp,'mwFreq')
  error('Exp.mwFreq is required for the thyme method.')
end

Opt.SinglePointDetection = false;

if predefinedExperiment
  % set up Exp structure from predefined experiment for thyme
  Exp = s_predefinedexperiments(Exp);
end


Pulses = cellfun(@isstruct,Exp.Sequence);
PulsePositions = find(Pulses);

% check if pulses have a finite length
idealPulses = false;

% and wether their definition requires Exp.mwFreq to be given
mwFreqrequired = false;

msgtp = 'Please provide pulse lenghts (Par.tp) for the following pulses:';
msgFrequency = 'Please provide Par.Frequency for the following pulses:';
for Pos = PulsePositions
  PulseNumber = find(PulsePositions==Pos);
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
  errmsg = ['Real pulses are required for the thyme-method.\n' msgtp '\n' msgFrequency];
elseif idealPulses
  errmsg = ['Real pulses are required for the thyme-method.\n' msgtp];
elseif mwFreqrequired
  errmsg = ['Real pulses are required for the thyme-method.\n' msgFrequency];
else
  errmsg = [];
end

if ~isempty(errmsg)
  error(sprintf(errmsg));
end


% Set up detection
if isfield(Exp,'DetWindow')
  % Validate Exp.DetWindow
  logmsg(1,'  setting up detection window:');
  if isfield(Exp,'DetSequence')
    warning('You provided Exp.DetWindow and Exp.DetSequence. Exp.DetSequence will be ignored.')
  end
  
  % Ensure that detection window does not overlap with a pulse
  if any(Exp.DetWindow<0) && isstruct(Exp.Sequence{end})
    error('You provided Exp.DetWindow with a negative value, but the last element in Exp.Sequence is a pulse. Detection into a pulse is not possible. Please adapt your detection or append an appropriate free evolution event.')
  elseif min(Exp.DetWindow) < 0 && (abs(min(Exp.DetWindow)) > Exp.Sequence{end})
    error('Your detection window is extending beyond the last free evolution into a pulse. Please shorten detection window or adapt the length of the free evolution.')
  end
  
  % Set up Exp.DetSequence - this requires adding a detection event of
  % Exp.Sequence
  if isstruct(Exp.Sequence{end})
    if Exp.DetWindow(1)>0
      Exp.Sequence{end+1} = Exp.DetWindow(1);
    end
  else
    Exp.Sequence{end} = Exp.Sequence{end} + Exp.DetWindow(1);
  end
  if length(Exp.DetWindow) == 1 || Exp.DetWindow(1) == Exp.DetWindow(2)
    % single point detection
    Opt.SinglePointDetection = true;
    Exp.Sequence{end+1} = 0;
    logmsg(1,'  single point detection');
  else
    % transient
    Exp.Sequence{end+1} = diff(Exp.DetWindow);
    logmsg(1,'  transient detection');
  end
  
  Exp.DetSequence = zeros(1,length(Exp.Sequence));
  Exp.DetSequence(end) = true;
  
elseif isfield(Exp,'DetSequence')
  % setting up detection in case of Exp.DetSequence
  logmsg(1,'  setting up detection:');
  if ischar(Exp.DetSequence) || isstring(Exp.DetSequence)
    % parsing strings
    if strcmp(Exp.DetSequence,'last')
      Exp.DetSequence = zeros(1,length(Exp.Sequence));
      Exp.DetSequence(end) = true;
      logmsg(1,'  detection is active during the last element in Exp.Sequence');
    elseif strcmp(Exp.DetSequence,'all')
      Exp.DetSequence = ones(1,length(Exp.Sequence));
      logmsg(1,'  all elements in Exp.Sequence are detected');
    else
      msg = 'The string you provided in Exp.DetSequence was not recognized';
      error(msg);
    end
  else
    if length(Exp.DetSequence) ~= 1 && length(Exp.DetSequence) ~= length(Exp.Sequence)
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
  Exp.DetSequence = ones(1,length(Exp.Sequence));
end

% Check if resonator is provided
IncludeResonator = false;
if any([isfield(Exp,'ResonatorFrequency') isfield(Exp,'ResonatorQL') isfield(Exp,'FrequencyResponse')])
  logmsg(1,'  resonator present');
  if isfield(Exp,'ResonatorFrequency') && isfield(Exp,'ResonatorQL')
    logmsg(1,['  resonator frequency: ' num2str(Exp.ResonatorFrequency)]);
    logmsg(1,['  resonator QL: ' num2str(Exp.ResonatorQL)]);
    Resonator.Arg1 = Exp.ResonatorFrequency;
    Resonator.Arg2 = Exp.ResonatorQL;
    IncludeResonator = true;
  elseif isfield(Exp,'ResonatorFrequency')
    error('Exp.ResonatorFrequency provided, but Exp.ResonatorQL is missing')
  elseif isfield(Exp,'ResonatorQL')
    error('Exp.ResonatorQL provided, but Exp.ResonatorFrequency is missing')
  end
  
  if isfield(Exp,'FrequencyResponse')
    if IncludeResonator
      logmsg(1,'  also found FrequencyResponse of the resonator, ignoring ResonatorFrequency and QL');
    end
    logmsg(1,'  using Exp.FrequencyResponse to simulate resonator.');
    Resonator.Arg1 = Exp.FrequencyResponse(1,:);
    Resonator.Arg2 = Exp.FrequencyResponse(2,:);
    IncludeResonator = true;
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

FreqShift = 0;

% get the absolute frequency shift, which is a combination of Exp.mwFreq
% (if given) and the shift of the simulation frame. this frequency shift is
% required for determining the minimum required step size.
if isfield(Exp,'mwFreq')
  FreqShift = FreqShift + Exp.mwFreq;
end

logmsg(1,'  setting up simulation frame:');
if isfield(Opt,'SimFreq')
  % user provided value fo shift of simulation frequency
  if Opt.SimFreq == 0
    % user requested lab frame
    FreqShift = FreqShift; %#ok<ASGSL> % Do nothing
    FrameShift = false;
  else
    FreqShift = FreqShift - Opt.SimFreq;
    FrameShift = Opt.SimFreq;
  end
  
else
  % determine minimum frequency in the experiment definition in order to
  % guess a frequency for the frame shift
  
  if isfield(Exp,'mwFreq') && Exp.mwFreq ~= 0
    minPulseFreq = Exp.mwFreq;
  else
    minPulseFreq = [];
  end
  
  % loop over the pulses and look for the minimum frequency there
  for iEvent = 1 : length(Exp.Sequence)
    if isstruct(Exp.Sequence{iEvent}) && isfield(Exp.Sequence{iEvent},'Frequency')
      if isempty(minPulseFreq)
        minPulseFreq = min(Exp.Sequence{iEvent}.Frequency/1000);
      else
        minPulseFreq = min([minPulseFreq (Exp.Sequence{iEvent}.Frequency/1000 + FreqShift)]);
      end
    end
  end
  
  % Have at least 2 GHz difference to the lowest frequency (for the frame
  % shift)
  FrameShift = floor(minPulseFreq-2);
  % only shift down, no upshifting - if lab frame frequencies exist that
  % are < 2 GHz
  if FrameShift > 0
    FreqShift = FreqShift - FrameShift;
  end
end

if FrameShift ~= 0
  logmsg(1,'  simulation frame frequency is %d GHz',FrameShift);
else
  logmsg(1,'  simulating in the lab frame');
end

Opt.FrameShift = FrameShift;
Opt.FreqShift = FreqShift;

% Check if IntTimeStep exists and if it is sufficient or, if none provided,
% compute a new one
logmsg(1,'  determining minimal required time step');
maxPulseFreq = FreqShift;
for iEvent = 1 : length(Exp.Sequence)
  if isstruct(Exp.Sequence{iEvent}) && isfield(Exp.Sequence{iEvent},'Frequency')
    maxPulseFreq = max([(Exp.Sequence{iEvent}.Frequency/1000 + FreqShift) maxPulseFreq]);
  end
end

Nyquist = 2*maxPulseFreq;
NyquistTime = 1/Nyquist/1000; % Time Step is in microseconds and Frequencies in GHz

IntStepFactor = 50; % control the size of the integration time step - 1/x of Nyquist
DetStepFactor = IntStepFactor/2; % gets the scaling factor for the detection time step (approx. 1/2 the size of the Nyquist step), in terms of multiples of integration time step
RefIntTimeStep = round(NyquistTime/IntStepFactor,2,'significant'); % The default integration time step
RefDetTimeStep = DetStepFactor*RefIntTimeStep; % The default dection time step

% validate time step
if isfield(Opt,'IntTimeStep') && ~isfield(Exp,'DetTimeStep')
  
  if Opt.IntTimeStep > NyquistTime
    errMsg = ['Your integration time step (Opt.IntTimeStep) does not fulfill the Nyquist criterion for the pulses you provided. For best results adapt it to ' num2str(RefIntTimeStep, '%10.1e') ' us or less.'];
    error(errMsg);
  elseif Opt.IntTimeStep > NyquistTime/(1/2*IntStepFactor)
    warnMsg = ['Although your integration time step (Opt.IntTimeStep) fulfills the Nyquist criterion for the pulses you provided, it might not be small enough for accurate results. You might want to adapt it to ' num2str(RefIntTimeStep, '%10.1e') ' us or less.'];
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
    errMsg = ['Your integration time step (Opt.IntTimeStep) does not fulfill the Nyquist criterion for the pulses you provided. For best results adapt it to ' num2str(RefIntTimeStep, '%10.1e') ' us or less.'];
    error(errMsg);
  elseif Opt.IntTimeStep > NyquistTime/(2/3*IntStepFactor)
    warnMsg = ['Although your integration time step (Opt.IntTimeStep) fulfills the Nyquist criterion for the pulses you provided, it might not be small enough for accurate results. You might want to adapt it to ' num2str(RefIntTimeStep, '%10.1e') ' us or less.'];
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

% Variables for bookkeeping of pulses and free evolution events
iPulse = 0;
iDelay = 0;

% A vector to quickly identify pulses, required for the reordering if
% pulses cross during the sequence
isPulse = zeros(1,length(Events));
PulseIndices = [];

Intervals = zeros(1,length(Exp.Sequence));

% making sure that relaxation is defined globally of for all individual
% elements in Exp.Sequence
if isfield(Opt,'Relaxation')
  if length(Opt.Relaxation) ~= 1 && length(Opt.Relaxation) ~= length(Exp.Sequence)
    error('The lengths of Exp.Sequence and Opt.Relaxation do not match. Length of Opt.Relaxation has to be 1 or the same as Exp.Sequence.')
  end
end

% making sure that state trajectories are defined globally of for all
% individual elements in Exp.Sequence
if isfield(Opt,'StateTrajectories')
  if length(Opt.StateTrajectories) ~= 1 && length(Opt.StateTrajectories) ~= length(Exp.Sequence)
    error('The lengths of Exp.Sequence and Opt.StateTrajectories do not match. Length of Opt.StateTrajectories has to be 1 or the same as Exp.Sequence.')
  end
end

% Setting up data structures for the pulses and events
for iEvent = 1 : nEvents
  if isstruct(Exp.Sequence{iEvent})
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

if isfield(Exp,'mwPolarization')
  if ~ischar(Exp.mwPolarization) && ~isstring(Exp.mwPolarization)
    error('Exp.mwPolarization has to be ''linear'' or ''circular''.')
  end
end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create the Eventstructure
% -------------------------------------------------------------------------
logmsg(1,'  parsing Exp.Sequence:');
if isfield(Exp,'DetWindow')
  logmsg(1,'  found %d pulse(s), %d free evolution period(s) and a detection window',nPulses,iDelay-1);
else
  logmsg(1,'  found %d pulse(s) and %d free evolution period(s)',nPulses,iDelay);
end
if nPulses > 0
  logmsg(1,'  computing wave forms and setting up the event structures');
else
  logmsg(1,'  setting up the event structures');
end

for iEvent = 1 : length(Exp.Sequence)
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
    
    if ~isfield(Exp.Sequence{iEvent},'IQ') % Pulse.IQ is being used for userdefined IQs
      % ---------------------------------------------------------------------
      % Pulse Specific Fields
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
        elseif length(Exp.Sequence{iEvent}.Frequency) > 1 && (Exp.Sequence{iEvent}.Frequency(2) ~= Exp.Sequence{iEvent}.Frequency(1))
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
      
      Pulse.PhaseCycle = ThisPhaseCycle;
      
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
        if IncludeResonator
          % if resonator is requested, pulses are elongated due to ringing.
          % the duration of ringing is stored in an additional field
          tOrig = t(end);
          [t,IQ] = resonator(t,IQ,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
          Events{iEvent}.Ringing = t(end) - tOrig;
        end
        % Shifts IQ of the pulse if necessary...
        if FreqShift ~= 0
          Opt.dt = Opt.IntTimeStep;
          [t, IQ] = rfmixer(t,IQ,FreqShift,'IQshift',Opt);
        end
        % ... and stores it in the event structure
        Events{iEvent}.IQ(iPCstep,:) = IQ;
      end
      
    else
      % in case user provided their own IQ(s)
      
      % determine input format:
      if iscell(Exp.Sequence{iEvent}.IQ)
        UserIQ = Exp.Sequence{iEvent}.IQ{1};
      elseif ismatrix(Exp.Sequence{iEvent}.IQ)
        UserIQ = Exp.Sequence{iEvent}.IQ;
      else
        errMsg = ['The data structure of the userdefined IQ of pulse on position ' num2str(iEvent) ' in Exp.Sequence is not recognized.'];
        error(errMsg);
      end
      
      [d1, d2] = size(UserIQ);
      % look for phase cycle and if found, verify that the IQ contains the
      % phase cycle
      if d1 ~= nPhaseSteps && d2 ~= nPhaseSteps
        errMsg = ['The dimensionality of the IQ signal provided of the pulse on position ' num2str(iEvent) ' in Exp.Sequence  is not in agreement with the phasecycle for this pulse. For user-defined waveforms the array must contain all IQs.'];
        error(errMsg);
      elseif d2 == nPhaseSteps
        UserIQ =  UserIQ';
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
        if IncludeResonator
          % but first if resonator is requested, pulses are elongated due to ringing.
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
      
      Pulse.userIQ.IQ = Exp.Sequence{iEvent}.IQ;
      Pulse.PhaseCycle = ThisPhaseCycle;
      
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
    if isfield(Exp,'mwPolarization') && ~isempty(Exp.mwPolarization) && strcmp(Exp.mwPolarization,'circular')
      Events{iEvent}.ComplexExcitation = true;
    else
      Events{iEvent}.ComplexExcitation = false;
    end
    
    % Temporarily store pulse paramaters to avoid reassigning them for creating the
    % vary table
    Pulse.EventIndex = iEvent;
    Pulses{iPulse} = Pulse;
    iPulse = iPulse + 1;
    
    Intervals(iEvent) = t(end);
  else
    % ---------------------------------------------------------------------
    % Delay/Free Evolution Specific Fields
    % ---------------------------------------------------------------------
    Events{iEvent}.type = 'free evolution';
    Events{iEvent}.t = Exp.Sequence{iEvent};
    
    Intervals(iEvent) = Exp.Sequence{iEvent};
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
  if ~isfield(Exp,'DetSequence') || isempty(Exp.DetSequence)
    Events{iEvent}.Detection = false;
  else
    if length(Exp.DetSequence) == 1
      Events{iEvent}.Detection = Exp.DetSequence;
    else
      Events{iEvent}.Detection = Exp.DetSequence(iEvent);
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
  Events{iEvent}.TimeStep = Opt.IntTimeStep;
  
  % Keep track of the frameshift
  Events{iEvent}.FrameShift = FrameShift;
  
end

% -------------------------------------------------------------------------
% Checks for overlap of pulses that are subject to ringing
% -------------------------------------------------------------------------
if IncludeResonator
  logmsg(1,'  checking for pulse overlap due to ringing from resonator:');
  
  for iEvent = PulseIndices
    FollowingEvent = iEvent + 1;
    if FollowingEvent <= length(Exp.Sequence) && strcmp(Events{FollowingEvent}.type,'pulse')
      error('When using a resonator, pulses need to be separated by inter pulse delays to accomodate for ringing from the resonator.')
    elseif FollowingEvent <= length(Exp.Sequence)
      ShortenedDelay = Events{FollowingEvent}.t - Events{iEvent}.Ringing;
      if ShortenedDelay < 0
        Msg = ['Event ' num2str(FollowingEvent) ' (a delay) is too short to accomodate ringing of the preceding pulse.'];
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
% structure can therefore be reduced
IncrementationScheme = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(Exp,'nPoints')
  logmsg(1,'-validating indirect dimensions------------------------');
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
  logmsg(1,'  found %d indirect dimension(s)',nDimensions);
  logmsg(1,'  with a total %d data points',prod(Vary.Points));
  for iDimension = 1 : nDimensions
    logmsg(1,'  parsing dimension no. %d',iDimension);
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
          if length(Strings) ~= 2 || ~strcmp(Strings{2},'Frequency') || length(Exp.(Field2Get){iLine,2}) == 1 || ~any(size(Exp.(Field2Get){iLine,2},1) == [1 Exp.nPoints(iDimension)-1])
            % determined if field was p1.Frequency
            if length(Exp.(Field2Get){iLine,2}) ~= 1 && length(Exp.(Field2Get){iLine,2}) ~= Exp.nPoints(iDimension)
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
                  PulseModifications{PulseNumber} = {iDimension Field Exp.(Field2Get){iLine,2} FieldIndex};
                else
                  n = size(PulseModifications{PulseNumber},1);
                  PulseModifications{PulseNumber}{n+1,1} = iDimension;
                  PulseModifications{PulseNumber}{n+1,2} = Field;
                  PulseModifications{PulseNumber}{n+1,3} = Exp.(Field2Get){iLine,2};
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
  logmsg(1,'  creating the Vary structure that contains all required wave forms and delay changes');
  for iDataPoint = 1 : nDataPoints
    % Load starting values for pulses and event lengths
    Pulses = InitialPulses;
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
          phase0 = Pulses{iPulse}.Phase;
          for iPCstep = 1 : size(Pulses{iPulse}.PhaseCycle,1)
            Pulses{iPulse}.Phase = phase0 + Pulses{iPulse}.PhaseCycle(iPCstep,1);
            [t,IQ] = pulse(Pulses{iPulse});
            if IncludeResonator
              % if a resonator is present, the ringing duration of each pulse
              % needs to be stored in the vary structure too
              tOrig = t(end);
              [t,IQ] = resonator(t,IQ,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
              Vary.Pulses{iPulse}.Ringing(ArrayIndex{:}) = t(end) - tOrig;
            end
            if FreqShift ~= 0
              Opt.dt = Opt.IntTimeStep;
              [~, IQ] = rfmixer(t,IQ,FreqShift,'IQshift',Opt);
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
          
          % get phase cycle
          nPhaseSteps = size(Pulses{iPulse}.PhaseCycle,1);
          [d1, d2] = size(UserIQ);
          
          % validate phase cycle
          if d1 ~= nPhaseSteps && d2 ~= nPhaseSteps
            errMsg = ['The dimensionality of the IQ signal provided for event ' num2str(iEvent) ' is not in agreement with the phasecycle for this pulse. For user-defined waveforms the array must contain all IQs.'];
            error(errMsg);
          elseif d2 == nPhaseSteps
            UserIQ =  UserIQ';
          end
          
          % Shifts IQ of the pulse if necessary...
          ShiftedUserIQ = [];
          for iPhaseStep = 1 : nPhaseSteps
            Opt.dt = Opt.IntTimeStep;
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

varargout{1} = Events;
varargout{2} = Vary;
varargout{3} = Opt;
varargout{4} = Exp;

logmsg(1,'  pulse sequence parsed successfully!');

end

function newstring = addtomessage(oldstring,toadd)
if isempty(oldstring)
  newstring = toadd;
else
  newstring = [oldstring ', ' toadd];
end
end