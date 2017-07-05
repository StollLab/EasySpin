function [Events, Vary, Opt] = sequencer(Exp,Opt)
% This function creates the Event and Vary Structures

Vary = [];

% Check if Timestep is sufficient
MaxFreq = max(abs(Exp.Frequency(:)));
Nyquist = 2*MaxFreq;

% Write this warning in a nicer way
if Exp.TimeStep > 1/Nyquist
  warning('Your Time Step (Exp.TimeStep) appears to not fullfill the Nyquist criterium for the provided frequencies.')
end

Events = cell(1,length(Exp.t));

iPulse = 1;
iDelay = iPulse;
DelayIndeces = zeros(1,length(Exp.t));
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
    
    % Looks for an excitation operator in the pulse definition structure
    % and if none is found, Sx is assumed
    if ~isfield(Pulse,'xOp')
      if isfield(Pulse,'ComplexExcitation') && Pulse.ComplexExcitation
        Events{iEvent}.xOp = spops(1/2,'x')+spops(1/2,'y');
      else
        Events{iEvent}.xOp = spops(1/2,'x');
      end
    else
      Events{iEvent}.xOp = Pulse.xOp;
    end
    
    % Checks if ComplexExcitation is requested for this Pulse, if not
    % specified Complex Excitation is switched off by default
    if ~isfield(Pulse,'ComplexExcitation')
      Events{iEvent}.ComplexExcitation = false;
    else
      Events{iEvent}.ComplexExcitation = Pulse.ComplexExcitation;
    end
    
    % Temporarily store pulse paramaters to avoid reassigning them for creating the
    % vary table
    Pulse.EventIndex = iEvent;
    Pulses{iPulse} = Pulse;
    
    % Incremeant the index for Pulse by 1
    iPulse = iPulse + 1;
    
  else
    % ---------------------------------------------------------------------
    % Delay/Free Evolution Specific Fields
    % ---------------------------------------------------------------------
    Events{iEvent}.type = 'free evolution';
    Events{iEvent}.t = 0:Exp.TimeStep:Exp.t(iEvent);
    DelayIndeces(iDelay) = iEvent;
    iDelay = iDelay + 1;
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
  if ~isfield(Opt,'Detection')
    if iEvent == length(Exp.t)
      Events{iEvent}.Detection = true;
    else
      Events{iEvent}.Detection = false;
    end
  else
    if length(Opt.Detection) == 1
      Events{iEvent}.Detection = Opt.Detection;
    elseif iEvent > length(Opt.Detection)
      error('You did not specify detection after Event %d.',iEvent);
    else 
      Events{iEvent}.Detection = Opt.Detection(iEvent);
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
  
  Events{iEvent}.Propagation = [];
    
end

if isfield(Exp,'nPoints')
  nDims = length(Exp.nPoints);
  Vary.Points = Exp.nPoints;
  
  for iDim = 1 : nDims
    
    if nDims == 1
      if isfield(Exp,'Inc')
        field2get = 'Inc';
      else
        field2get = 'Inc1';
      end
    else
      field2get = ['Inc' num2str(iDim)];
    end
    
    VariedEvents = zeros(1,size(Exp.(field2get),1));
    iModified = 1;
    for iLines = 1 : size(Exp.(field2get),1)
      FullString = Exp.(field2get){iLines,1};

      SplitStrings = regexp(FullString,',','split');
      
      for iModifiedEvent = 1 : length(SplitStrings)
              
        Strings = regexp(SplitStrings{iModifiedEvent},'\.','split');
        
        type = Strings{1}(1);
        index = str2double(Strings{1}(2:end));
        
        
        switch type
          case 'p'
            field = Strings{2};
            
            EventNumber = Pulses{index}.EventIndex;
            PulseNumber = index;
            
            VariedEvents(iModified) = EventNumber;
            Pulse = Pulses{PulseNumber};
            
            
            Vary.IQs{EventNumber}{1} = Events{EventNumber}.IQ;
            Vary.ts{EventNumber}{1} = Events{EventNumber}.t;
            
            switch field
              case 't'
                field = 'tp';
            end
            
            Start = Pulse.(field);
            
            for iPoint = 2 : Vary.Points(iDim)
              Pulse.(field) = Start + (iPoint-1)*Exp.(field2get){iLines,2};
              for iPCstep = 1 : size(Pulse.PhaseCycle,1)
                Pulse.Phase = Pulse.Phase+Pulse.PhaseCycle(iPCstep,1);
                [t,IQ] = pulse(Pulse);
                if iPCstep == 1
                  IQs = zeros(size(Pulse.PhaseCycle,1),length(IQ));
                end
                IQs(iPCstep,:) = IQ;
              end
              
              Vary.IQs{EventNumber}{iPoint} = IQs;
              Vary.ts{EventNumber}{iPoint} = t;
            end
            
          case 'd'
            
            EventNumber = DelayIndeces(index);
            
            VariedEvents(iModified) = EventNumber;
            
            Start = Exp.t(EventNumber);
            Vary.ts{EventNumber}{1} = Events{EventNumber}.t;
            for iPoint = 2 : Vary.Points(iDim)
              t = Start + (iPoint-1)*Exp.(field2get){iLines,2};
              
              Vary.ts{EventNumber}{iPoint} = 0:Exp.TimeStep:t;
            end
            
        end
        iModified = iModified + 1;
      end
      
      
    end
    Vary.Events{iDim} = VariedEvents;
  end
else
  %make empty vary table here
end







