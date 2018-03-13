% spidyan    Simulate spindyanmics during pulse EPR experiments
%
%     [t,Signal] = spidyan(Sys,Exp,Opt)
%     [t,Signal,out] = spidyan(Sys,Exp,Opt)
%
%     Inputs:
%       Sys   ... spin system with electron spin and ESEEM nuclei
%       Exp   ... experimental parameters (time unit us)
%       Opt   ... simulation options
%
%     Outputs:
%       t                 ... time axis
%       Signal            ... simulated signals of detected events
%       out               ... output structure with fields:
%         .FinalState        ... density matrix/matrices at the end of the
%                                experiment
%         .StateTrajectories ... cell array with density matrices from each
%                                timestep during evolution
%         .Events            ... structure containing events

function varargout = spidyan(Sys,Exp,Opt) 

if (nargin==0), help(mfilename); return; end

% Input argument scanning, get display level and prompt
%=======================================================================
% Check Matlab version
VersionErrorStr = chkmlver;
error(VersionErrorStr);

% --------License ------------------------------------------------
LicErr = 'Could not determine license.';
Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------

% Guard against wrong number of input or output arguments.
if (nargin<2) || (nargin>3), error('Wrong number of input arguments!'); end
if (nargout<0), error('Not enough output arguments.'); end
if (nargout>3), error('Too many output arguments.'); end

% Initialize options structure to zero if not given.
if (nargin<3), Opt = struct('unused',NaN); end
if isempty(Opt), Opt = struct('unused',NaN); end

if ~isstruct(Sys)
  error('First input argument (Sys) must be a structure!');
end
if ~isstruct(Exp)
  error('Second input argument (Exp) must be a structure!');
end
if ~isstruct(Opt)
  error('Third input argument (Opt) must be a structure!');
end
% 

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;


%----------------------------------------------------------------------
% Preprocess Input
%----------------------------------------------------------------------

% Check for magnetic field
if ~isfield(Exp,'Field')
  error('Exp.Field is required.')
elseif length(Exp.Field) ~= 1
  error('Exp.Field must be a single number (in mT).');
end

% Adapt Zeeman frequencies and g values for the selected simulation frame, 
% if they are provided
if isfield(Opt,'FrameShift') && ~isempty(Opt.FrameShift)
  if Opt.FrameShift < 0
    warning('The value provided for Opt.FrameShift is negative, but should be positive. The simulation frequency is always shifted to lower frequency and hence does not require a negative sign.')
    Opt.FrameShift = abs(Opt.FrameShift);
  end
  if isfield(Sys,'ZeemanFreq')
    Sys.ZeemanFreq =  Sys.ZeemanFreq - Opt.FrameShift;
  end
  
  if isfield(Sys,'g')  
    nElectrons = length(Sys.S);
    
    gshift = (Opt.FrameShift*1e9)*planck/bmagn/(Exp.Field(end)*1e-3);
    
    issize = @(A,siz) all(size(A)==siz);
    fullg = issize(Sys.g,[3*nElectrons 3]);
    if fullg
      gshiftmat = repmat(gshift*eye(3),[nElectrons,1]);
      Sys.g = Sys.g - gshiftmat;
    else
      Sys.g = Sys.g - gshift;
    end
    
  end
end

Exp.Frequency = Exp.Frequency*1000; % GHz -> MHz

% Translate Frequency to g values
if isfield(Sys,'ZeemanFreq')
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

% Remove field ZeemanFreq if given, which is spidyan specific
if isfield(Sys,'ZeemanFreq')
  Sys = rmfield(Sys,'ZeemanFreq');
end

% Build the Event and Vary structure
[Events, Vary] = s_sequencer(Exp,Opt);


% Validate and build spin system as well as excitation operators
[Sys, Sigma, DetOps, Events, Relaxation] = s_propagationsetup(Sys,Events,Opt);

% Get Hamiltonian
Ham = sham(Sys,Exp.Field*[0 0 1]);

%----------------------------------------------------------------------
% Propagation
%----------------------------------------------------------------------

% Calls the actual propagation engine
[TimeAxis, RawSignal, FinalState, StateTrajectories, Events] = ...
  s_thyme(Sigma, Ham, DetOps, Events, Relaxation, Vary);

%----------------------------------------------------------------------
% Signal Processing
%----------------------------------------------------------------------

% Signal postprocessing, such as down conversion and filtering and
% checking output of the timeaxis 
if ~isempty(RawSignal)
  
  if ~iscell(TimeAxis)
    if size(unique(TimeAxis,'rows'),1) == 1
      TimeAxis = TimeAxis(1,:);
    end
  else
    TimeAxisChanged = false;
    for iCell = 2 : length(TimeAxis)
      if ~isequal(TimeAxis{1},TimeAxis{iCell})
        TimeAxisChanged = true;
        break
      end
    end
    
    if ~TimeAxisChanged
      TimeAxis = TimeAxis{1};
    end
    
  end

  
  % Adapt FreqTranslation if needed
  if ~isfield(Opt,'FreqTranslation') || isempty(Opt.FreqTranslation)
    FreqTranslation = [];
    
  else
    nDetOps = numel(DetOps);
        
    if isfield(Opt,'FrameShift') && ~isempty(Opt.FrameShift)
      Opt.FreqTranslation(Opt.FreqTranslation > 0) = Opt.FreqTranslation(Opt.FreqTranslation > 0) - Opt.FrameShift;
      Opt.FreqTranslation(Opt.FreqTranslation < 0) = Opt.FreqTranslation(Opt.FreqTranslation < 0) + Opt.FrameShift;
    end
    FreqTranslation = zeros(1,nDetOps);
    FreqTranslation(1:length(Opt.FreqTranslation)) = Opt.FreqTranslation;
    
  end
  
  % Downconversion/processing of signal
  Signal = signalprocessing(TimeAxis,RawSignal,FreqTranslation);

else
  Signal = [];
end

if size(FinalState,1) == 1
  FinalState = squeeze(FinalState);
end

% Assigning outputs
switch nargout
  case 1
    varargout = {Signal};
  case 2
    varargout = {TimeAxis,Signal};
  case 3
    out.FinalState = FinalState;
    out.StateTrajectories = StateTrajectories;
    out.Events = Events;
    varargout = {TimeAxis,Signal,out};
  otherwise
    error('Incorrect number of output arguments. 1,2, or 3 expected.');
end
