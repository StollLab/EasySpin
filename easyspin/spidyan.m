% spidyan    Simulate spindyanmics during pulse EPR experiments
%
%     [t,Signal] = spidyan(Sys,Exp,Opt)
%     [t,Signal,info] = spidyan(Sys,Exp,Opt)
%
%     Inputs:
%       Sys   ... spin system with electron spin and nuclei
%       Exp   ... experimental parameters (time unit µs)
%       Opt   ... simulation options
%
%     Outputs:
%       t                 ... time axis
%       Signal            ... simulated signals of detected events
%       info              ... structure with fields:
%         .FinalState        ... density matrix/matrices at the end of the
%                                experiment
%         .StateTrajectories ... cell array with density matrices from each
%                                timestep during evolution
%         .Events            ... structure containing events

function varargout = spidyan(Sys,Exp,Opt) 

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version
warning(chkmlver);

StartTime = clock;

% Input argument scanning, get display level and prompt
%-------------------------------------------------------------------------------

% Guard against wrong number of input or output arguments.
if nargin<2 || nargin>3, error('Wrong number of input arguments!'); end
if nargout<0, error('Not enough output arguments.'); end
if nargout>3, error('Too many output arguments.'); end

% Initialize options structure to zero if not given.
if nargin<3, Opt = struct; end
if isempty(Opt), Opt = struct; end

if ~isstruct(Sys)
  error('First input argument (Sys) must be a structure!');
end
if ~isstruct(Exp)
  error('Second input argument (Exp) must be a structure!');
end
if ~isstruct(Opt)
  error('Third input argument (Opt) must be a structure!');
end

% Set log level
if isfield(Opt,'Verbosity'), eslogger(Opt.Verbosity); end

eslogger(1,'=begin=spidyan====%s=================',char(datetime));
eslogger(2,'  log level %d',eslogger);

% Preprocess Input
%-------------------------------------------------------------------------------

% Check for magnetic field
if ~isfield(Exp,'Field')
  eslogger(1,'  no Exp.Field given');
  eslogger(1,'  checking if Sys allows for running without a field...');
  if isfield(Sys,'g')
    error('Exp.Field is required for Sys.g.')
  elseif Sys.ZeemanFreq % && (~isfield(Sys,'S') || all(Sys.S==1/2))
    % if Sys.ZeemanFreq is used, it Exp.Field is not required, however the
    % program needs a reference field to correctly shift frequencies and
    % to calculate nuclear Larmor frequencies (if nuclei are given). The
    % program does so, by "guessing" a field from Sys.ZeemanFreq and gfree
    Exp.Field = (Sys.ZeemanFreq(1)*1e6)*planck/bmagn/(gfree*1e-6);
    eslogger(1,'  Exp.Field is needed for calculating simulation frame frequencies');
    eslogger(1,'  setting Exp.Field to %.1f mT',Exp.Field);
  else
    error('Exp.Field is required for given spin system.')
  end  
elseif length(Exp.Field) ~= 1
  error('Exp.Field must be a single number (in mT).');
end

% Build the Event and Vary structure
Opt.SimulationMode = 'thyme'; % to force s_sequencer to setup for step wise

% If the user provides relaxation times, but no Opt.Relaxation,
% relaxationis switched on for the entire sequence
if (isfield(Sys,'T1') || isfield(Sys,'T2')) && ~isfield(Opt,'Relaxation')
  Opt.Relaxation = true;
end

[Events, Vary, Opt] = s_sequencer(Exp,Opt);

% check if output can be plotted in case nargout == 0
if nargout == 0 && isfield(Exp,'nPoints')
  if Opt.SinglePointDetection && length(Exp.nPoints) > 2
    error('spidyan can not plot more than two indirect dimensions in combination with single point detection.');
  elseif ~Opt.SinglePointDetection && length(Exp.nPoints) > 1
    error('spidyan can not plot more than one indirect dimension in combination with transient detection.');
  end
end

% Adapt Zeeman frequencies for the selected simulation frame
eslogger(1,'-validating spin system--------------------------------');

if ~isfield(Sys,'S')
 if isfield(Sys,'ZeemanFreq') && ~isscalar(Sys.ZeemanFreq)
  error('If you provide more than one value in Sys.ZeemanFreq, Sys.S must be specified');
 end
end

if isfield(Sys,'ZeemanFreq') && Opt.FrameShift ~= 0
  eslogger(1,'adapting Sys.ZeemanFreq to the simulation frame');
  Sys.ZeemanFreq =  Sys.ZeemanFreq - Opt.FrameShift;
end

% Adapt g-values for the selected simulation frame
if isfield(Sys,'g')
  if isfield(Sys,'S')
    nElectrons = length(Sys.S);
  else
    nElectrons = 1;
  end
  
  if Opt.FrameShift~=0
    eslogger(1,'  adapting Sys.g to the simulation frame');
  end
  gshift = (Opt.FrameShift*1e9)*planck/bmagn/(Exp.Field(end)*1e-3);
  
  fullg = isequal(size(Sys.g),[3*nElectrons 3]);
  if fullg
    gshiftmat = repmat(gshift*eye(3),[nElectrons,1]);
    Sys.g = Sys.g - gshiftmat;
  else
    Sys.g = Sys.g - gshift;
  end
end

% Translate Frequency to g values
if isfield(Sys,'ZeemanFreq')
  eslogger(1,'  translating Sys.ZeemanFreq into g values');
  Sys.ZeemanFreq = Sys.ZeemanFreq*1000; % GHz -> MHz
  if isfield(Sys,'g')
    [~, dgTensor] = size(Sys.g);
  else
    dgTensor = 1; 
  end
  % Recalculates the g value
  for ieSpin = 1 : length(Sys.ZeemanFreq)
    if Sys.ZeemanFreq(ieSpin) ~= 0
      g = (Sys.ZeemanFreq(ieSpin)*1e6)*planck/bmagn/(Exp.Field*1e-3);
      Sys.g(ieSpin,1:dgTensor) = g;
    end
  end
end

% Remove spidyan-specific field ZeemanFreq if given
if isfield(Sys,'ZeemanFreq')
  eslogger(1,'  removing Sys.ZeemanFreq from Sys structure');
  Sys = rmfield(Sys,'ZeemanFreq');
end

% Transfer DetOperator from Exp to Opt, to make it available to
% s_propagationsetup
if isfield(Exp,'DetOperator')
  % catch the case when Exp.DetOperator = 'z1' instead of {'z1'}
  if ~iscell(Exp.DetOperator)
    Exp.DetOperator = {Exp.DetOperator};
  end
  Opt.DetOperator = Exp.DetOperator;
end

if isfield(Exp,'DetPhase')
  if isfield(Opt,'DetOperator') && length(Exp.DetPhase) ~= length(Opt.DetOperator)
    error('Exp.DetPhase has to have the same length as Exp.DetOperator.');
  elseif ~isfield(Opt,'DetOperator') && length(Exp.DetPhase) ~= 1
    error('Exp.DetPhase has to have length one if no Exp.DetOperator is provided.');
  end  
end
  

% Validate and build spin system as well as excitation operators
eslogger(1,'  parsing the Sys structure...');

[Sys, Sigma, DetOps, Events, Relaxation] = s_propagationsetup(Sys,Events,Opt);

% Get Hamiltonian
eslogger(1,'  computing lab frame Hamiltonian');
H = ham(Sys,Exp.Field*[0 0 1]);

% Propagation
%-------------------------------------------------------------------------------
if ~isempty(Relaxation)
  eslogger(1,'  adapting relaxation superoperator to system frame');
  [U,~] = eig(H);
  R = kron(transpose(U),(U'));
  Relaxation.Gamma = R'*Relaxation.Gamma*R;
end

% Calls the actual propagation engine
eslogger(1,'-starting propagation----------------------------------');
eslogger(1,'  depending on the Exp setup this may take a while...');
[TimeAxis, RawSignal, FinalState, StateTrajectories, Events] = ...
  s_thyme(Sigma, H, DetOps, Events, Relaxation, Vary);
eslogger(1,'  propagation finished!');

% Signal Processing
%-------------------------------------------------------------------------------
eslogger(1,'-validating and processing outout----------------------');
% Signal postprocessing, such as down conversion and filtering and
% checking output of the timeaxis 

nDetOps = numel(DetOps); % number of detection operators

% Applying detection phase
if ~isempty(RawSignal) && isfield(Exp,'DetPhase')
  eslogger(1,'  applying detection phase: %d*pi',Exp.DetPhase/pi);
  DetPhase = exp(-1i*Exp.DetPhase);
  
  RawSignal = applyPhaseShift(RawSignal,DetPhase);
end

if Opt.SinglePointDetection
  eslogger(1,'  single point detection...');
  if isfield(Exp,'nPoints') && nDetOps ~= 1
    DimSignal =  ndims(RawSignal);
    Signal = permute(RawSignal,[1:(DimSignal-1) DimSignal+1 DimSignal]);
  else
    Signal = RawSignal;
  end
else
  if ~isempty(RawSignal)
    % Check if Exp.DetFreq and Exp.DetOperator have the same length
    if isfield(Exp,'DetFreq')
      eslogger(1,'  veryfing Exp.DetFreq...');
      if isfield(Exp,'DetOperator')
        if length(Exp.DetFreq) ~= length(Exp.DetOperator)
          error('Exp.DetOperator and Exp.DetFreq must contain the same number of elements.')
        end
      elseif length(Exp.DetFreq) ~= 1
        error('If Exp.DetOperator is not provided the default detection operator (S+ for all electrons) is assumed. In this case Exp.DetFreq must contain only one frequency.')
      end
    else
      if ~isfield(Exp,'DetOperator') && isfield(Exp,'mwFreq')
        eslogger(1,'  using Exp.mwFreq as detection frequency');
        Exp.DetFreq = Exp.mwFreq;
      end
    end
    
    eslogger(1,'  processing transients...');
    
    % Adapt FreqTranslation if needed
    FreqTranslation = zeros(1,nDetOps);
    
    if isfield(Exp,'DetFreq') && ~isempty(Exp.DetFreq)
      
      % This adapts the values for DetFreq to simulation frame
      Exp.DetFreq(Exp.DetFreq > 0) = Exp.DetFreq(Exp.DetFreq > 0) - Events{1}.FrameShift;
      Exp.DetFreq(Exp.DetFreq < 0) = Exp.DetFreq(Exp.DetFreq < 0) + Events{1}.FrameShift;
      
      % And then writes them
      FreqTranslation(1:length(Exp.DetFreq)) = - Exp.DetFreq; % To make it a down conversion for negative frequencies add the neg sign
      
    end
    
    % Downconversion/processing of signal
    Signal = signalprocessing(TimeAxis,RawSignal,FreqTranslation);
    
    % If time axis is the same for each data point, it is reduced to a single
    % vector at this point - helps with plotting
    if ~iscell(TimeAxis)
      SizeT = size(TimeAxis);
      linearTimeAxis = reshape(TimeAxis,[prod(SizeT(1:end-1)) SizeT(end)]);
      if size(unique(linearTimeAxis,'rows'),1) == 1
        TimeAxis = linearTimeAxis(1,:);
      end
    end
  else
    eslogger(1,'  nothing was detected...');
    Signal = [];
  end
end

% Reduce the dimensionality of the final state for simulations with only
% one acquisition point
if ndims(FinalState) == 3 && size(FinalState,1) == 1
  FinalState = squeeze(FinalState);
end

% Same for the StateTrajectories (if any exist). If only one vector of
% StateTrajectories was recorded, the nested structure is removed
if ~isempty(StateTrajectories)
  SizeStateTrajectories = size(StateTrajectories);
  if all(SizeStateTrajectories == 1)
    StateTrajectories = StateTrajectories{1};
  end
end

% Assigning outputs
switch nargout
  case 0
    % no output argument - graphical output
    eslogger(1,'-no output requested------------------------------------');
    eslogger(1,'  switching to graphical output');
    if isempty(Signal)
      disp('Detection was switched off, nothing to display.')
    else
      s_plotting(TimeAxis,Signal,Exp,Opt);
    end
  case 1
    varargout = {Signal};
  case 2
    varargout = {TimeAxis,Signal};
  case 3
    info.FinalState = FinalState;
    info.StateTrajectories = StateTrajectories;
    info.Events = Events;
    varargout = {TimeAxis,Signal,info};
  otherwise
    error('Incorrect number of output arguments. 1,2, or 3 expected.');
end

% Report performance
%-------------------------------------------------------------------------------
hmsString = elapsedtime(StartTime,clock);
eslogger(1,'spidyan took %s',hmsString);

eslogger(1,'=end=spidyan======%s=================\n',char(datetime));

end

%===============================================================================
function phasedSignal = applyPhaseShift(rawSignal, phase)

nPhases = length(phase);

if iscell(rawSignal)
  nPoints = numel(rawSignal);
  
  phasedSignal = cell(size(rawSignal));
  
  for iPoint = 1:nPoints
    phasedSignal{iPoint} = bsxfun(@times,rawSignal{iPoint},phase');
  end
  
else
  signalSize = size(rawSignal);
  
  nPoints = prod(signalSize(1:end-2));
  % reshape the n-dimensional array into a 3-dimensional array, that
  % can be looped over linearly along the first dimension (which
  % corresponds to all the acquistion points)
  rawSignal = reshape(rawSignal,[nPoints signalSize(end-1) signalSize(end)]);
  
  phasedSignal = zeros(size(rawSignal));
  
  for iPhase = 1 : nPhases
    phasedSignal(:,iPhase,:) = rawSignal(:,iPhase,:)*phase(iPhase);
  end
  
  phasedSignal = reshape(phasedSignal,signalSize);
  
end

end
