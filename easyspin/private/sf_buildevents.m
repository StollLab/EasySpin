% sf_buildevents  Build events and indirect dimensions for the full engine
%
%   [Events,Vary,Opt] = sf_buildevents(Exp,Opt)
%
%   Parses the pulse sequence specified in Exp, using settings from Opt, to
%   create the Events and Vary structures needed for the simulation by
%   sf_propagate (full engine of saffron). Predefined experiments
%   are expanded by saffron before.
%
%   Input:
%     Exp  experiment structure, with fields
%       Sequence            pulse sequence: cell array of pulses and delays (µs)
%       DetWindow           detection window (one or two numbers, in µs,
%                           relative to the end of the sequence), or
%                           'sequence' for the whole sequence
%       nPoints             number of points of the indirect dimensions
%       Dim1, Dim2, ...     incrementation scheme
%       mwFreq              microwave carrier frequency (GHz), optional
%       rfFreq              RF carrier frequency (MHz), optional
%       DetTimeStep         detection time step (µs), optional
%       mwPolarization      'linear' (default) or 'circular'
%       ResonatorFrequency, ResonatorQL, FrequencyResponse, ResonatorMode
%                           resonator, see resonator()
%       pulseCycles         phases of each pulse in the phase cycle, cell
%                           array with one element per pulse, one phase
%                           (radians, added to Pulse.Phase) per row;
%                           optional, default 0
%     Opt  options structure, with fields
%       rfGammaRef          gyromagnetic ratio (MHz/mT) for nominal RF flip
%                           angles
%       rfEnhancement       true (default) to include the hyperfine
%                           enhancement of RF pulses
%       maxTransitionFreq   upper bound on the transition frequencies (MHz)
%                           of the spin systems; needed for RF pulses in the
%                           lab frame with rfEnhancement (see rfTimeStep)
%       SimFreq             simulation frame frequency (GHz), optional
%       IntTimeStep         integration time step (µs), optional
%       Relaxation          true/false, or one value per element of
%                           Exp.Sequence, to turn on relaxation during
%                           the events
%
%   Output:
%     Events  cell array of structures, one per event (the elements of
%             Exp.Sequence, plus events appended for the detection window),
%             with the fields
%       type                'pulse' or 'free evolution'
%       t                   pulses: time axis (µs); free evolution: length (µs)
%       Relaxation          true if relaxation is on during this event
%       Detection           true if the signal is detected during this event
%       StateTrajectories   store the density matrices during the event
%                           (always false here; set by callers of
%                           sf_propagate, e.g. tests)
%       Propagation         empty, filled in by sf_propagate
%       TimeStep            time step (µs) of the integration
%       DetTimeStep         detection time step (µs)
%       FrameShift          simulation frame frequency (GHz)
%      and, for pulses only,
%       IQ                  waveform in the simulation frame, one row per
%                           phase cycle step; ideal pulses: Flip*exp(1i*phase)
%       Ideal               true for ideal pulses
%       Channel             'mw' or 'rf'
%       PhaseCycle          phases (radians), one row per phase cycle step
%                           (from Exp.pulseCycles)
%       FrameOffset         carrier frequency (GHz) in the simulation frame;
%                           for RF pulses, see rfFrameParameters
%       FrameFreq, FrameBand  RF pulses only: in a shifted simulation frame,
%                           frequency (MHz) of the nuclear RF frame, and
%                           [lowest highest] RF frequency (MHz) of the pulse;
%                           empty in the lab frame
%       Polarization        'linear' or 'circular'
%       Ringing             duration (µs) of the resonator ringing, if any
%     Vary    structure with the indirect dimensions, [] if there are none,
%             with the fields
%       Points              number of points of each dimension
%       IncrementationTable cell array with one element per dimension: array
%                           of the changes in event lengths (µs), one row
%                           per event and one column per point; empty if
%                           no event lengths change
%       Pulses              cell array with one element per pulse, empty if
%                           the pulse is not varied; otherwise, a structure
%                           with cell arrays over the points of IQs and ts
%                           (waveforms and time axes), array Ringing (if
%                           any), and, for RF pulses, FrameOffsets,
%                           FrameFreqs, FrameBands and TimeSteps
%     Opt     input structure, with the fields
%       SinglePointDetection  true for single-point detection
%       FrameShift            simulation frame frequency (GHz)
%       IntTimeStep           integration time step (µs)


function varargout = sf_buildevents(Exp,Opt)

logmsg(1,'-parsing pulse sequence-----------------------------');

Vary = [];

% Without Exp.mwFreq, the microwave carrier is at 0: pulse frequencies
% (Pulse.Frequency, in MHz) are absolute
noCarrier = ~isfield(Exp,'mwFreq');
% The same for RF pulses: carrier Exp.rfFreq (MHz), or 0 without it (checked
% by saffron)
if isfield(Exp,'rfFreq')
  rfCarrier = Exp.rfFreq;
else
  rfCarrier = 0;
end

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
% below (its length is checked by saffron)
nUserEvents = numel(Exp.Sequence);

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
  % Update Exp.Sequence based on Exp.DetWindow (saffron checks that the
  % window does not overlap with a pulse)
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
    
    % Phases of this pulse in the phase cycle, one per row
    if isfield(Exp,'pulseCycles') && ~isempty(Exp.pulseCycles)
      thisPhaseCycle = Exp.pulseCycles{iPulse};
    else
      thisPhaseCycle = 0;
    end
    nPhaseSteps = size(thisPhaseCycle,1);
    
    if ~isfield(Exp.Sequence{iEvent},'IQ') % Pulse.IQ is being used for userdefined IQs
      % ---------------------------------------------------------------------
      % Pulse specific fields
      % ---------------------------------------------------------------------
      
      % if no pulse type is provided, the default is a rectangular pulse
      if ~isfield(Exp.Sequence{iEvent},'Type')
        
        % First check for field frequency and correct it
        if ~isfield(Exp.Sequence{iEvent},'Frequency')
          % if no frequency is defined, frequency is set to 0 (the carrier
          % will be added later)
          Exp.Sequence{iEvent}.Frequency = 0;
          Exp.Sequence{iEvent}.Type = 'rectangular';
          if ~isfield(Exp,'mwFreq') && ~isRFEvent(iEvent)
            error('If you do not provide a frequency for the pulse, you need to give Exp.mwFreq.')
          end
        elseif ~isscalar(Exp.Sequence{iEvent}.Frequency) && (Exp.Sequence{iEvent}.Frequency(2) ~= Exp.Sequence{iEvent}.Frequency(1))
          % make sure that, if a frequency was provided, it is not a frequency sweep
          msg = ['Pulse at position ' num2str(iEvent) ' in Exp.Sequence: no Pulse.Type specified, assuming a monochromatice rectangular pulse, but the field Pulse.Frequency looks like a frequency-swept pulse.'];
          error(msg);
        else
          Exp.Sequence{iEvent}.Type = 'rectangular';
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
      elseif ~isfield(Pulse,'Qcrit') && ~isfield(Pulse,'Amplitude') && ~isfield(Pulse,'B2')
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
        Pulse.rfCarrier = rfCarrier;  % MHz
        Pulse.rfFrame = frameShift~=0;
        Pulse.Detected = detectEvents(iEvent);
        % In the lab frame with the hyperfine enhancement, the RF operator
        % includes the electron magnetic moment, so the sampling of the RF
        % field must not drive the electron transitions (see rfTimeStep)
        labEnhancement = ~Pulse.rfFrame && (~isfield(Opt,'rfEnhancement') || Opt.rfEnhancement);
        if labEnhancement
          Pulse.maxTransitionFreq = Opt.maxTransitionFreq;  % MHz
        else
          Pulse.maxTransitionFreq = 0;
        end
        Pulse.TimeStep = rfTimeStep(Pulse,Exp.DetTimeStep);
        if ~isfield(Pulse,'B2') && isfield(Opt,'rfGammaRef')
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
      % User-defined waveform (Pulse.IQ): a single waveform, or a cell
      % array with one waveform per point of an indirect dimension (the
      % first one is used here)
      if iscell(Exp.Sequence{iEvent}.IQ)
        userIQ = Exp.Sequence{iEvent}.IQ{1};
      else
        userIQ = Exp.Sequence{iEvent}.IQ;
      end
      checkUserIQ(userIQ,iEvent);

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
      
      % Wave forms of all steps of the phase cycle
      [t,Events{iEvent}.IQ,ringing] = userWaveforms(tIQ,userIQ,thisPhaseCycle,Opt,includeResonator,Resonator,Exp);
      if ~isempty(ringing)
        Events{iEvent}.Ringing = ringing;
      end
      
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
    % RF pulses have their own carrier (Exp.rfFreq, or 0 without it, in
    % GHz), on the same sequence clock. In a shifted simulation frame, they
    % are propagated in a nuclear frame rotating at FrameFreq (see
    % rfFrameParameters).
    if Pulse.isRF
      [Events{iEvent}.FrameOffset,Events{iEvent}.FrameFreq,Events{iEvent}.FrameBand] = ...
        rfFrameParameters(Pulse);
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
    % Lookup table of the propagators (see pulsePropagators in sf_propagate)
    % with a finer amplitude resolution than for microwave pulses: RF
    % pulses stepped with the lab-frame field have thousands of steps, and
    % the rounding of the field to 1024 levels gave errors of up to 2e-4
    % relative to the echo, which depended on the RF phase.
    Events{iEvent}.vertRes = 8192;
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
            pulses{iPulse}.TimeStep = rfTimeStep(pulses{iPulse},Exp.DetTimeStep);
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
            [Vary.Pulses{iPulse}.FrameOffsets{ArrayIndex{:}},...
              Vary.Pulses{iPulse}.FrameFreqs{ArrayIndex{:}},...
              Vary.Pulses{iPulse}.FrameBands{ArrayIndex{:}}] = rfFrameParameters(pulses{iPulse});
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
            IndexToSave = ones(1,length(DimensionIndices));
            IndexToSave(Dimension) = DimensionIndices(Dimension);
            ArrayIndex = num2cell(IndexToSave);
            
            userIQ = pulses{iPulse}.userIQ.IQ{IndexToLoad};
            
            tIQ = pulses{iPulse}.userIQ.t{IndexToLoad};
          else
            IndexToLoad = ones(1,length(DimensionIndices));
            
            for iModification = 1 : size(PulseModifications{iPulse},1)
              Dimension = PulseModifications{iPulse}{iModification,1};
              IndexToLoad(Dimension) = DimensionIndices(Dimension);
            end
            ArrayIndex = num2cell(IndexToLoad);
            
            userIQ = pulses{iPulse}.userIQ.IQ{ArrayIndex{:}};
            
            tIQ = pulses{iPulse}.userIQ.t{ArrayIndex{:}};
            
          end
          
          % Wave forms of all steps of the phase cycle
          checkUserIQ(userIQ,pulses{iPulse}.EventIndex);
          [t,Vary.Pulses{iPulse}.IQs{ArrayIndex{:}},ringing] = userWaveforms(tIQ,userIQ,...
            pulses{iPulse}.PhaseCycle,Opt,includeResonator,Resonator,Exp);
          if ~isempty(ringing)
            Vary.Pulses{iPulse}.Ringing(ArrayIndex{:}) = ringing;
          end
          
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
% pulseWaveform  Wave form of a pulse for one step of its phase cycle
%
%   [t,IQ,ringing] = pulseWaveform(Pulse,iPCstep,Opt,freqShift,includeResonator,Resonator,Exp)
%
%   Computes the wave form in the simulation frame:
%   - Ideal microwave pulses: t = 0 and a single complex number
%     IQ = Flip*exp(1i*phase), from which sf_propagate builds the
%     instantaneous rotation.
%   - RF pulses: lab-frame waveform (carrier at Pulse.rfCarrier, MHz, from
%     the start of the pulse, plus the pulse's Frequency, which is absolute
%     if the carrier is 0) in units of the RF field (mT). In a shifted
%     simulation frame (Pulse.rfFrame), the envelope relative to the center
%     frequency f0 of Pulse.Frequency instead, i.e. without the carrier and
%     with exp(-2i*pi*f0*t). Without B2, the amplitude is nominal: the flip
%     angle of a bare nucleus with the gyromagnetic ratio Pulse.rfGammaRef
%     (MHz/mT).
%   - Finite microwave pulses: from pulse(), through the resonator, shifted
%     to the simulation frame.
%
%   Input:
%     Pulse             pulse parameters (see pulse()), with additional fields
%       isIdeal, isRF     pulse type
%       PhaseCycle        phases (radians), one row per phase cycle step
%       rfCarrier         RF carrier frequency (MHz, 0 without a carrier),
%                         RF pulses
%       rfFrame           true for the envelope in the nuclear RF frame,
%                         RF pulses
%       B2                RF field amplitude (mT), RF pulses, optional
%       rfGammaRef        gyromagnetic ratio (MHz/mT), RF pulses without B2
%     iPCstep           index of the phase cycle step
%     Opt               options; uses IntTimeStep
%     freqShift         frequency (GHz) of the microwave carrier in the
%                       simulation frame
%     includeResonator  true to include the resonator
%     Resonator         resonator parameters, fields Arg1, Arg2, Arg3 are
%                       the arguments of resonator()
%     Exp               experiment; uses mwFreq
%
%   Output:
%     t        time axis (µs); [0 0] for ideal pulses
%     IQ       wave form; single number for ideal pulses
%     ringing  duration (µs) of the resonator ringing, empty without resonator
function [t,IQ,ringing] = pulseWaveform(Pulse,iPCstep,Opt,freqShift,includeResonator,Resonator,Exp)
ringing = [];
phase = Pulse.Phase + Pulse.PhaseCycle(iPCstep);
if Pulse.isIdeal
  t = [0 0];  % zero duration (two time points for detection bookkeeping)
  IQ = Pulse.Flip*exp(1i*phase);
  return
end
Par = Pulse;
Par.Phase = phase;
Par = rmfield(Par,intersect(fieldnames(Par),{'isIdeal','isRF','rfCarrier','rfGammaRef','maxTransitionFreq',...
  'rfFrame','Detected','B2','Channel','PhaseCycle','EventIndex','ArrayIndex','userIQ','Operator'}));
if Pulse.isRF
  if isfield(Pulse,'B2') && ~isempty(Pulse.B2)
    Par = rmfield(Par,intersect(fieldnames(Par),{'Flip','Qcrit'}));
    Par.Amplitude = Pulse.B2;  % mT
    [t,IQ] = pulse(Par);
  else
    if ~isfield(Pulse,'rfGammaRef')
      error('RF pulses in the full engine need B2 (in mT), or a nucleus to define a nominal flip angle.');
    end
    [t,IQ] = pulse(Par);
    IQ = IQ/Pulse.rfGammaRef;  % MHz -> mT
  end
  if Pulse.rfFrame
    % envelope relative to the center frequency f0 of the pulse
    IQ = IQ.*exp(-2i*pi*mean(Pulse.Frequency([1 end]))*t);
  else
    IQ = IQ.*exp(2i*pi*Pulse.rfCarrier*t);
  end
  return
end
[t,IQ] = pulse(Par);
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
% userWaveforms  Wave forms of a user-defined pulse for all phase cycle steps
%
%   [t,IQ,ringing] = userWaveforms(tIQ,userIQ,PhaseCycle,Opt,includeResonator,Resonator,Exp)
%
%   Each step of the phase cycle shifts the phase of the waveform. The
%   waveform, upconverted with the carrier at Exp.mwFreq, is passed
%   through the resonator and shifted to the simulation frame.
%
%   Input:
%     tIQ               time axis of the waveform (µs)
%     userIQ            waveform (Pulse.IQ), vector
%     PhaseCycle        phases (radians), one row per step
%     Opt               options; uses IntTimeStep, FrameShift
%     includeResonator  true to include the resonator
%     Resonator         resonator parameters, fields Arg1, Arg2, Arg3 are
%                       the arguments of resonator()
%     Exp               experiment; uses mwFreq
%
%   Output:
%     t        time axis (µs)
%     IQ       wave forms in the simulation frame, one row per step
%     ringing  duration (µs) of the resonator ringing, empty without resonator
function [t,IQ,ringing] = userWaveforms(tIQ,userIQ,PhaseCycle,Opt,includeResonator,Resonator,Exp)
ringing = [];
Opt.dt = Opt.IntTimeStep;
nSteps = size(PhaseCycle,1);
IQ = [];
for iStep = 1:nSteps
  tStep = tIQ;
  IQstep = exp(1i*PhaseCycle(iStep))*userIQ(:).';
  if includeResonator
    % the pulse is elongated by the ringing
    [tStep,IQstep] = resonator(tStep,IQstep,Exp.mwFreq,Resonator.Arg1,Resonator.Arg2,Resonator.Arg3);
    ringing = tStep(end) - tIQ(end);
  end
  [t,IQstep] = rfmixer(tStep,IQstep,-Opt.FrameShift,'IQshift',Opt);
  if iStep==1, IQ = zeros(nSteps,numel(IQstep)); end
  IQ(iStep,:) = IQstep;
end
end

%-------------------------------------------------------------------------------
% checkUserIQ  Check a user-defined waveform
%
%   checkUserIQ(userIQ,iEvent)
%
%   Input:
%     userIQ  waveform (Pulse.IQ)
%     iEvent  position of the pulse in Exp.Sequence
function checkUserIQ(userIQ,iEvent)
if ~isnumeric(userIQ) || ~isvector(userIQ)
  error(['Pulse.IQ of the pulse at position %d in Exp.Sequence must be a single waveform (a vector), ' ...
    'or a cell array of waveforms for an indirect dimension. To change its phase in a phase cycle, ' ...
    'use Exp.PhaseCycle.'],iEvent);
end
end

%-------------------------------------------------------------------------------
% rfTimeStep  Time step for RF pulses
%
%   dt = rfTimeStep(Pulse,DetTimeStep)
%
%   Lab-frame waveform: at least 40 points per period of the highest RF
%   frequency, and at least 100 steps per pulse (needed for RF frequencies
%   near zero). With the electron magnetic moment in the RF operator
%   (hyperfine enhancement), the field, constant within each step, also
%   has components at k/dt +- f_RF (k = 1,2,...; f_RF the highest RF
%   frequency), with amplitudes of about f_RF*dt/k relative to the RF
%   field. Where they meet an electron transition, they drive it (spurious
%   ENDOR lines). So the lowest of them, 1/dt - f_RF, is kept above all
%   transition frequencies, with a margin of 20%:
%   1/dt >= 1.2*(Pulse.maxTransitionFreq + f_RF). The propagator of each
%   step is exact, so the electron precession itself needs no finer steps.
%   Envelope in the nuclear RF frame (Pulse.rfFrame): the field is constant
%   within each step, so the step only needs to resolve the envelope. A
%   rectangular pulse at a fixed frequency is a single step. Otherwise, at
%   least 40 points per period of the largest offset of Pulse.Frequency
%   from its center, and at least 100 steps per pulse (amplitude shapes).
%   The step is at most DetTimeStep if the pulse is detected.
%   The step divides the pulse length, so that the waveform ends at tp.
%
%   Input:
%     Pulse        pulse; uses tp, Frequency (offsets from the carrier, or
%                  absolute frequencies without a carrier, MHz), Type,
%                  rfCarrier (MHz, 0 without a carrier), rfFrame, Detected,
%                  maxTransitionFreq (MHz, upper bound on the transition
%                  frequencies of the spin systems, 0 without hyperfine
%                  enhancement)
%     DetTimeStep  detection time step (µs)
%
%   Output:
%     dt           time step (µs)
function dt = rfTimeStep(Pulse,DetTimeStep)
tp = Pulse.tp;
if ~Pulse.rfFrame
  fRF = abs(Pulse.rfCarrier) + max(abs(Pulse.Frequency));  % MHz
  nSteps = max(100,ceil(40*fRF*tp));
  nSteps = max(nSteps,ceil(1.2*(Pulse.maxTransitionFreq + fRF)*tp));
  dt = tp/nSteps;
  return
end
maxOffset = max(abs(Pulse.Frequency - mean(Pulse.Frequency([1 end]))));
amplitudeShape = strtok(Pulse.Type,'/');
if strcmpi(amplitudeShape,'rectangular') && maxOffset==0
  nSteps = 1;
else
  nSteps = max(100,ceil(40*maxOffset*tp));
end
if Pulse.Detected
  nSteps = max(nSteps,ceil(tp/DetTimeStep*(1-1e-12)));
end
dt = tp/nSteps;
end

%-------------------------------------------------------------------------------
% rfFrameParameters  Carrier phase and nuclear frame of an RF pulse
%
%   [FrameOffset,FrameFreq,FrameBand] = rfFrameParameters(Pulse)
%
%   The waveform gets the start phase 2*pi*FrameOffset*1e3*tStart, tStart
%   being the start time of the pulse on the sequence clock.
%   - Lab-frame waveform: FrameOffset is the carrier, so that the carrier
%     phase is continuous over the sequence. sf_propagate multiplies the
%     waveform with exp(1i*startPhase).
%   - Envelope in the nuclear RF frame (Pulse.rfFrame), relative to the
%     center frequency f0 of Pulse.Frequency: the lab field is the envelope
%     times exp(2i*pi*FrameFreq*t) on the sequence clock, with FrameFreq =
%     rfCarrier + f0. With FrameOffset = -f0, the pulse waveform starts at
%     the start of the pulse, as in the lab frame. sf_propagate takes the
%     field relative to the start of the pulse and adds the frame phase
%     2*pi*FrameFreq*tStart to the start phase; rfFramePropagators applies
%     the start phase as a rotation of the whole pulse.
%
%   Input:
%     Pulse        pulse; uses Frequency (MHz), rfCarrier (MHz), rfFrame
%
%   Output:
%     FrameOffset  frequency (GHz) of the start-time phase
%     FrameFreq    frequency (MHz) of the nuclear RF frame; empty in the
%                  lab frame
%     FrameBand    [lowest highest] RF frequency of the pulse (MHz); empty
%                  in the lab frame
function [FrameOffset,FrameFreq,FrameBand] = rfFrameParameters(Pulse)
if ~Pulse.rfFrame
  FrameOffset = Pulse.rfCarrier/1e3;
  FrameFreq = [];
  FrameBand = [];
  return
end
f0 = mean(Pulse.Frequency([1 end]));
FrameOffset = -f0/1e3;
FrameFreq = Pulse.rfCarrier + f0;
FrameBand = Pulse.rfCarrier + [min(Pulse.Frequency) max(Pulse.Frequency)];
end
