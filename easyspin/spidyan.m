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

% Build the Event and Vary structure
[Events, Vary, Opt] = s_sequencer(Exp,Opt);

% Adapt Zeeman frequencies and g values for the selected simulation frame, 
% if they are provided
if isfield(Sys,'ZeemanFreq')
  Sys.ZeemanFreq =  Sys.ZeemanFreq - Opt.FrameShift;
end

if isfield(Sys,'g')
  if isfield(Sys,'S')
    nElectrons = length(Sys.S);
  else
    nElectrons = 1;
  end
  
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

% Validate and build spin system as well as excitation operators
if isfield(Exp,'DetOperator')
  Opt.DetOperator = Exp.DetOperator;
end

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
if Opt.SinglePointDetection
  if isfield(Exp,'nPoints')
    DimSignal =  ndims(RawSignal);
    Signal = permute(RawSignal,[1:(DimSignal-1) DimSignal+1 DimSignal]);
  else
    Signal = RawSignal;
  end
else
  if ~isempty(RawSignal)
    
    % Adapt FreqTranslation if needed
    nDetOps = numel(DetOps);
    FreqTranslation = zeros(1,nDetOps);
    
    if isfield(Exp,'DetFrequency') && ~isempty(Exp.DetFrequency)
      
      % This adapts the values for DetFrequency to simulation frame
      Exp.DetFrequency(Exp.DetFrequency > 0) = Exp.DetFrequency(Exp.DetFrequency > 0) - Events{1}.FrameShift;
      Exp.DetFrequency(Exp.DetFrequency < 0) = Exp.DetFrequency(Exp.DetFrequency < 0) + Events{1}.FrameShift;
      
      % And then writes them
      FreqTranslation(1:length(Exp.DetFrequency)) = - Exp.DetFrequency; % To make it a down conversion for negative frequencies add the neg sign
      
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
