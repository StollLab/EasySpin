% This function parses the pulse sequence specified in the Exp structure and
% additionally using settings from the Opt structure to create the Event and
% Vary structures needed for the simulation by sf_propagate (full-space
% engine of saffron). Predefined experiments (macros) are expanded by
% saffron before.
%
% Input:
%   Exp - structure with the following fields
%     Sequence
%     DetWindow  - detection window (one or two numbers, in µs, relative
%                  to the end of the sequence), or 'sequence' for the
%                  whole sequence
%     Dim1, Dim2, etc
%     mwFreq
%     DetFreq
%     DetTimeStep
%     nPoints
%     mwPolarization
%     ResonatorFrequency
%     ResonatorQL
%     FrequencyResponse
%     ResonatorMode
%     PhaseCycle
%   Opt - structure with the following fields
%     rfGammaRef - gyromagnetic ratio (MHz/mT) for nominal RF flip angles
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
%     Polarization
%     Relaxation
%     Detection
%     StateTrajectories - store the density matrices during the event
%         (always false here; set by callers of sf_propagate, e.g. tests)
%     Propagation
%     FrameShift
%     TimeStep
%     DetTimeStep
%   Vary - structure with the following fields
%     Points
%     IncrementationTable
%     Pulses
%   Opt - strucure with the following fields
%     SinglePointDetection
%     FrameShift
%     %FreqShift


function varargout = sf_buildevents(Exp,Opt)

logmsg(1,'-parsing pulse sequence-----------------------------');

Vary = [];

% Without Exp.mwFreq, the microwave carrier is at 0: pulse frequencies
% (Pulse.Frequency, in MHz) are absolute
noCarrier = ~isfield(Exp,'mwFreq');

Opt.SinglePointDetection = false;

pulsePositions = find(cellfun(@isstruct,Exp.Sequence));

% Pulse types: ideal microwave pulses (no tp and no IQ), RF pulses
% (Channel 'rf'), and finite microwave pulses
isRFEvent = false(1,numel(Exp.Sequence));
isIdealEvent = false(1,numel(Exp.Sequence));
for Pos = pulsePositions
  p_ = Exp.Sequence{Pos};
  if numel(p_)>1
    error('Simultaneous pulses (struct arrays in Exp.Sequence) are not supported yet.');
  end
  if isfield(p_,'Channel') && strcmpi(p_.Channel,'rf')
    isRFEvent(Pos) = true;
    if ~isfield(Exp,'rfFreq')
      error('Exp.rfFreq (RF carrier frequency, in MHz) is required for RF pulses.');
    end
    if ~isfield(p_,'tp') || p_.tp<=0
      error('RF pulses need a pulse length (tp>0).');
    end
  elseif ~isfield(p_,'IQ') && (~isfield(p_,'tp') || p_.tp==0)
    isIdealEvent(Pos) = true;
    if ~isfield(p_,'Flip')
      error('Ideal pulses (no tp) need a flip angle (Flip).');
    end
  elseif noCarrier && ~isfield(p_,'IQ') && ...
      (~isfield(p_,'Frequency') || isempty(p_.Frequency) || all(p_.Frequency==0))
    % Pulse.Frequency = 0 would be a pulse at 0 MHz in the lab frame
    error(['Without Exp.mwFreq, the frequencies of microwave pulses are absolute: ' ...
      'give a nonzero Pulse.Frequency (in MHz) for every finite pulse, or give Exp.mwFreq.']);
  end
  if isfield(p_,'IQ') && ~isfield(Opt,'IntTimeStep')
    error('If you use userprovided IQs for pulses, Opt.IntTimeStep must be given.')
  end
end
mwPulsePositions = pulsePositions(~isRFEvent(pulsePositions));
finitemwPulsePositions = pulsePositions(~isRFEvent(pulsePositions) & ~isIdealEvent(pulsePositions));

% Event-specific relaxation (Opt.Relaxation) refers to the elements of
% Exp.Sequence as given by the user, before detection events are appended
% below
nUserEvents = numel(Exp.Sequence);
if isfield(Opt,'Relaxation') && ~isscalar(Opt.Relaxation) && numel(Opt.Relaxation)~=nUserEvents
  error('Opt.Relaxation must be a single value or have one entry per element of Exp.Sequence (%d).',nUserEvents);
end

% Set up detection
%-------------------------------------------------------------------------------
logmsg(1,'  setting up detection');

% detectEvents: detection on (true) or off for each event
detectWholeSequence = ~isfield(Exp,'DetWindow') || ischar(Exp.DetWindow) || isstring(Exp.DetWindow);
if detectWholeSequence
  % Exp.DetWindow = 'sequence' (checked by saffron)
  logmsg(1,'  detection during the whole sequence');
  detectEvents = true(1,numel(Exp.Sequence));
  % a sequence of a single zero-length delay is a single point
  Opt.SinglePointDetection = isscalar(Exp.Sequence) && ~isstruct(Exp.Sequence{1}) && Exp.Sequence{1}==0;
else
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
  
  % Only the appended detection event is detected
  detectEvents = false(1,numel(Exp.Sequence));
  detectEvents(end) = true;
end

% Detection events appended for Exp.DetWindow are free evolution, so they
% take the relaxation setting of the last delay in the user's sequence (or,
% without delays, relaxation is on if it is on anywhere)
userDelays = find(~cellfun(@isstruct,Exp.Sequence(1:nUserEvents)));
if isfield(Opt,'Relaxation') && ~isscalar(Opt.Relaxation)
  if isempty(userDelays)
    flag = any(Opt.Relaxation);
  else
    flag = Opt.Relaxation(userDelays(end));
  end
  Opt.Relaxation(nUserEvents+1:numel(Exp.Sequence)) = flag;
end

% Check if resonator is provided
includeResonator = false;
Resonator = struct;
if any([isfield(Exp,'ResonatorFrequency') isfield(Exp,'ResonatorQL') isfield(Exp,'FrequencyResponse')])
  logmsg(1,'  resonator present');
  if noCarrier
    error('A resonator (Exp.ResonatorFrequency, Exp.FrequencyResponse) needs the microwave carrier frequency Exp.mwFreq.');
  end
  if isfield(Exp,'ResonatorFrequency') && isfield(Exp,'ResonatorQL')
    logmsg(1,'  resonator frequency: %s',num2str(Exp.ResonatorFrequency));
    logmsg(1,'  resonator QL: %s',num2str(Exp.ResonatorQL));
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
      logmsg(1,'  resonator mode: %s',Exp.ResonatorMode);
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
  % User-provided value for simulation frequency shift; frames rotate in
  % the positive sense (as Exp.DetFreq)
  if ~isscalar(Opt.SimFreq) || ~isreal(Opt.SimFreq) || Opt.SimFreq<0
    error('Opt.SimFreq must be a single nonnegative frequency (in GHz).');
  end
  frameShift = Opt.SimFreq;
  freqShift = freqShift - frameShift;
  
else

  % Determine minimum frequency used by microwave pulses (RF pulses don't
  % determine the frame)
  minFreq = inf;
  for iEvent = mwPulsePositions
    if isfield(Exp.Sequence{iEvent},'Frequency')
      pulseFreq = Exp.Sequence{iEvent}.Frequency/1e3;  % MHz -> GHz
    else
      pulseFreq = 0;
    end
    minFreq = min(minFreq,freqShift+min(pulseFreq));
  end
  
  % Have at least 2 GHz difference to the lowest frequency (for the frame shift)
  % Only shift down, not up - if lab frame frequencies exist that are < 2 GHz,
  % simulate in the lab frame
  limFreq = 2;  % GHz
  if isinf(minFreq), minFreq = freqShift; end  % no microwave pulses
  frameShift = max(floor(minFreq-limFreq),0);
  freqShift = freqShift - frameShift;
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

% Determine maximum frequency used by microwave pulses (relative to sim freq)
maxFreq = -inf;
for iEvent = mwPulsePositions
  if isfield(Exp.Sequence{iEvent},'Frequency')
    pulseFreq = Exp.Sequence{iEvent}.Frequency/1e3;  % MHz -> GHz
  else
    pulseFreq = 0;
  end
  maxFreq = max(max(abs(freqShift+pulseFreq)),maxFreq);  % pulses can be below the simulation frame
end

if isinf(maxFreq), maxFreq = abs(freqShift); end  % no microwave pulses
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
  
  % largest multiple of the integration time step that does not exceed
  % the default detection time step
  logmsg(1,'  automatically assuming a suitable detection time step');
  Exp.DetTimeStep = max(1,floor(RefDetTimeStep/Opt.IntTimeStep))*Opt.IntTimeStep;

elseif ~isfield(Opt,'IntTimeStep') && isfield(Exp,'DetTimeStep')
  
  logmsg(1,'  computing an integration time step that fits the provided Exp.DetTimeStep');
  factor = ceil(Exp.DetTimeStep/RefIntTimeStep);
  Opt.IntTimeStep = Exp.DetTimeStep/factor; 
  
elseif isfield(Opt,'IntTimeStep') && isfield(Exp,'DetTimeStep')
  
  ratio = Exp.DetTimeStep/Opt.IntTimeStep;
  if round(ratio)<1 || abs(ratio-round(ratio)) > 1e-9*ratio
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

if ~detectWholeSequence
  logmsg(1,'  found %d pulse(s), %d free evolution period(s) and a detection window',nPulses,nDelays-1);
else
  logmsg(1,'  found %d pulse(s) and %d free evolution period(s)',nPulses,nDelays);
end

if isfield(Exp,'mwPolarization') && ~isempty(Exp.mwPolarization)
  if ~any(strcmp(Exp.mwPolarization,{'linear','circular'}))
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
        Pulse.tp = 0;  % ideal pulse
      end
      if ~isfield(Pulse,'Frequency')
        Pulse.Frequency = [0 0];
      end

      Pulse.PhaseCycle = thisPhaseCycle;

      % Gets the flip angle
      if isfield(Exp.Sequence{iEvent},'Flip')
        Pulse.Flip = Exp.Sequence{iEvent}.Flip;
      elseif ~isfield(Pulse,'Qcrit') && ~isfield(Pulse,'Amplitude') && ~isfield(Pulse,'B1')
        error('Flip angle for pulse no. %d is missing.',iPulse)
      end

      % Gets the phase for the pulse, if none is provided, the phase is
      % assumed to be 0
      if ~isfield(Pulse,'Phase')
        Pulse.Phase = 0;
      end

      % Channel and time step
      Pulse.isIdeal = isIdealEvent(iEvent);
      Pulse.isRF = isRFEvent(iEvent);
      if Pulse.isRF
        Pulse.rfCarrier = Exp.rfFreq;  % MHz
        Pulse.TimeStep = rfTimeStep(Exp.rfFreq,Pulse.Frequency,Pulse.tp);
        if ~isfield(Pulse,'B1') && isfield(Opt,'rfGammaRef')
          Pulse.rfGammaRef = Opt.rfGammaRef;
        end
      else
        Pulse.TimeStep = Opt.IntTimeStep;
      end

      % Wave forms of all steps of the phase cycle
      for iPCstep = 1 : nPhaseSteps
        [t,IQ,ringing] = pulseWaveform(Pulse,iPCstep,Opt,freqShift,includeResonator,Resonator,Exp);
        if ~isempty(ringing)
          Events{iEvent}.Ringing = ringing;
        end
        Events{iEvent}.IQ(iPCstep,:) = IQ;
      end
      Events{iEvent}.Ideal = Pulse.isIdeal;
      Events{iEvent}.Channel = 'mw';
      if Pulse.isRF
        Events{iEvent}.Channel = 'rf';
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
      Pulse.isIdeal = false;
      Pulse.isRF = false;
      Events{iEvent}.Ideal = false;
      Events{iEvent}.Channel = 'mw';

    end
    
    % Store the time axis of the pulse in the Events structure
    Events{iEvent}.t = t;
    
    % Store the PhaseCycle in the Event structure
    Events{iEvent}.PhaseCycle = thisPhaseCycle;
    
    % Frequency (GHz) of the microwave carrier in the simulation frame. The
    % pulse shape was shifted into the simulation frame on the time axis of
    % the pulse, and sf_propagate adds the carrier phase at the start time of the
    % pulse, so that the carrier is continuous over the sequence. This
    % includes user-supplied lab-frame waveforms, which are taken to be
    % upconverted with the carrier at Exp.mwFreq.
    % RF pulses have their own carrier (absolute RF frequency, GHz), on the
    % same sequence clock.
    if Pulse.isRF
      Events{iEvent}.FrameOffset = Exp.rfFreq/1e3;
    else
      Events{iEvent}.FrameOffset = freqShift;
    end

    % Polarization of the microwave field (linear by default), used by
    % sf_propagationsetup to build the pulse Hamiltonian
    if Pulse.isRF
      Events{iEvent}.Polarization = 'linear';
    elseif isfield(Exp,'mwPolarization') && ~isempty(Exp.mwPolarization)
      Events{iEvent}.Polarization = char(Exp.mwPolarization);
    else
      Events{iEvent}.Polarization = 'linear';
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
  
  % Detection during this event (set up from Exp.DetWindow above)
  Events{iEvent}.Detection = detectEvents(iEvent);

  % Density matrices during the event are not stored
  Events{iEvent}.StateTrajectories = false;

  % Store an empty propagation structure, filled in by sf_propagate
  Events{iEvent}.Propagation = [];
  
  % Store the time steps, which will be needed in sf_propagate to calculate time
  % axis and propagators: pulses are integrated and detected with
  % TimeStep, detected free evolution is computed at steps of DetTimeStep
  Events{iEvent}.TimeStep = Opt.IntTimeStep;
  if isfield(Events{iEvent},'Channel') && strcmp(Events{iEvent}.Channel,'rf')
    Events{iEvent}.TimeStep = Pulse.TimeStep;  % RF waveform step
  end
  Events{iEvent}.DetTimeStep = Exp.DetTimeStep;
  
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
        % RF carrier frequency: changes the carrier of all RF pulses
        if strcmp(strtrim(SplitStrings{iModifiedEvent}),'rfFreq')
          rfPulses = find(cellfun(@(P)P.isRF,pulses));
          if isempty(rfPulses)
            error('Exp.%s changes rfFreq, but the sequence has no RF pulses.',dimX);
          end
          for iP = rfPulses
            newMod = {iDimension 'rfCarrier' Exp.(dimX){iLine,2} []};
            if isempty(PulseModifications{iP})
              PulseModifications{iP} = newMod;
            else
              PulseModifications{iP}(end+1,:) = newMod;
            end
          end
          continue
        end
        if any(strcmp(strtrim(SplitStrings{iModifiedEvent}),{'mwFreq','Field'}))
          error('Indirect dimensions over Exp.%s are not supported yet.',strtrim(SplitStrings{iModifiedEvent}));
        end
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

          % RF pulses: time step for the RF frequency of this point
          if pulses{iPulse}.isRF
            pulses{iPulse}.TimeStep = rfTimeStep(pulses{iPulse}.rfCarrier,pulses{iPulse}.Frequency,pulses{iPulse}.tp);
          end

          % Convert array into cell for indexing
          ArrayIndex = num2cell(pulses{iPulse}.ArrayIndex);
              
          % Compute wave form and store it
          for iPCstep = 1 : size(pulses{iPulse}.PhaseCycle,1)
            [t,IQ,ringing] = pulseWaveform(pulses{iPulse},iPCstep,Opt,freqShift,includeResonator,Resonator,Exp);
            if ~isempty(ringing)
              Vary.Pulses{iPulse}.Ringing(ArrayIndex{:}) = ringing;
            end
            % ... and stores it in the vary structure
            Vary.Pulses{iPulse}.IQs{ArrayIndex{:}}(iPCstep,:) = IQ;
          end

          Vary.Pulses{iPulse}.ts{ArrayIndex{:}} = t;
          if pulses{iPulse}.isRF
            Vary.Pulses{iPulse}.FrameOffsets{ArrayIndex{:}} = pulses{iPulse}.rfCarrier/1e3;
            Vary.Pulses{iPulse}.TimeSteps{ArrayIndex{:}} = pulses{iPulse}.TimeStep;
          end
          
          
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
    [newSequence, newEventLengths] = sf_reorder_events(EventLengths,isPulse);
    
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

varargout = {Events,Vary,Opt};

logmsg(1,'  pulse sequence parsed successfully!');

end

%-------------------------------------------------------------------------------
% Wave form of a pulse for one step of its phase cycle, in the simulation
% frame. Phase cycle rows are [phase weight] or [phase weight scale].
% - Ideal microwave pulses: t = 0 and a single complex number
%   IQ = Flip*exp(1i*phase), from which sf_propagate builds the instantaneous
%   rotation.
% - RF pulses: lab-frame waveform (carrier at Pulse.rfCarrier, MHz, from
%   the start of the pulse) in units of the RF field (mT). Without B1, the
%   amplitude is nominal: the flip angle of a bare nucleus with the
%   gyromagnetic ratio Pulse.rfGammaRef (MHz/mT).
% - Finite microwave pulses: from pulse(), through the resonator, shifted
%   to the simulation frame.
function [t,IQ,ringing] = pulseWaveform(Pulse,iPCstep,Opt,freqShift,includeResonator,Resonator,Exp)
ringing = [];
pc = Pulse.PhaseCycle;
if size(pc,2)<2, pc = [0 1]; end
phase = Pulse.Phase + pc(iPCstep,1);
scale = 1;
if size(pc,2)>=3, scale = pc(iPCstep,3); end
if Pulse.isIdeal
  t = [0 0];  % zero duration (two time points for detection bookkeeping)
  IQ = scale*Pulse.Flip*exp(1i*phase);
  return
end
Par = Pulse;
Par.Phase = phase;
Par = rmfield(Par,intersect(fieldnames(Par),{'isIdeal','isRF','rfCarrier','rfGammaRef',...
  'B1','Channel','PhaseCycle','EventIndex','ArrayIndex','userIQ','Operator'}));
if Pulse.isRF
  if isfield(Pulse,'B1') && ~isempty(Pulse.B1)
    Par = rmfield(Par,intersect(fieldnames(Par),{'Flip','Qcrit'}));
    Par.Amplitude = Pulse.B1;  % mT
    [t,IQ] = pulse(Par);
  else
    if ~isfield(Pulse,'rfGammaRef')
      error('RF pulses in the full-space engine need B1 (in mT), or a nucleus to define a nominal flip angle.');
    end
    [t,IQ] = pulse(Par);
    IQ = IQ/Pulse.rfGammaRef;  % MHz -> mT
  end
  IQ = scale*IQ.*exp(2i*pi*Pulse.rfCarrier*t);
  return
end
[t,IQ] = pulse(Par);
IQ = scale*IQ;
if includeResonator
  % pulses are elongated due to ringing; its duration is returned
  tOrig = t(end);
  [t,IQ] = resonator(t,IQ,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
  ringing = t(end) - tOrig;
end
if freqShift~=0
  Opt.dt = Opt.IntTimeStep;
  [t,IQ] = rfmixer(t,IQ,freqShift,'IQshift',Opt);
end
end

%-------------------------------------------------------------------------------
% Time step for RF pulses (µs): at least 40 points per period of the
% highest RF frequency in the lab frame, and at least 100 steps per pulse
% (needed for RF frequencies near zero). The step divides the pulse length,
% so that the waveform ends at tp.
function dt = rfTimeStep(rfFreq,Frequency,tp)
dt = tp/max(100,ceil(40*(abs(rfFreq)+max(abs(Frequency)))*tp));
end
