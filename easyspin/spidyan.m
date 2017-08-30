% spidyan    Simulate pulse EPR spectra
%
%     [x,S] = saffron(Sys,Exp,Opt)
%     [x,S,out] = saffron(Sys,Exp,Opt)
%
%     [x1,x2,S] = saffron(Sys,Exp,Opt)
%     [x1,x2,S,out] = saffron(Sys,Exp,Opt)
%
%     Sys   ... spin system with electron spin and ESEEM nuclei
%     Exp   ... experimental parameters (time unit us)
%     Opt   ... simulation options
%
%     out:
%       x       ... time or frequency axis (1D experiments)
%       x1, x2  ... time or frequency axis (2D experiments)
%       S       ... simulated signal (ESEEM) or spectrum (ENDOR)
%       out     ... structure with FFT of ESEEM signal

function [TimeAxis, Signal,FinalState,StateTrajectories,NewEvents] = spidyan(Sys,Exp,Det,Opt)

% if (nargin==0), help(mfilename); return; end

% Get time for performance report at the end.
% StartTime = clock;
% 
% % Input argument scanning, get display level and prompt
% %=======================================================================
% % Check Matlab version
% VersionErrorStr = chkmlver;
% error(VersionErrorStr);
% 
% % --------License ------------------------------------------------
% LicErr = 'Could not determine license.';
% Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------

% Guard against wrong number of input or output arguments.
% if (nargin<2) || (nargin>3), error('Wrong number of input arguments!'); end
% if (nargout<0), error('Not enough output arguments.'); end
% if (nargout>4), error('Too many output arguments.'); end

% Initialize options structure to zero if not given.
% if (nargin<3), Opt = struct('unused',NaN); end
% if isempty(Opt), Opt = struct('unused',NaN); end


% if ~isstruct(Exp)
%   error('Second input argument (Exp) must be a structure!');
% end
% if ~isstruct(Opt)
%   error('Third input argument (Opt) must be a structure!');
% end
% 

% A global variable sets the level of log display. The global variable
% is used in logmsg(), which does the log display.
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 0; end
global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;



if length(Exp.Field) == 1
  B = [0 0 Exp.Field];
else
  B = Exp.Field;
end

% Translate Frequency to g values
if isfield(Sys,'ZeemanFreq')
  if isfield(Sys,'g')
    [~, dgTensor] = size(Sys.g);
  else
    dgTensor = 1;
  end
  % Recalculates the g value 
  for ieSpin = 1 : length(Sys.ZeemanFreq)
    if Sys.ZeemanFreq(ieSpin) ~= 0
      g = Sys.ZeemanFreq(ieSpin)*1e9*planck/bmagn/Exp.Field(end);
      Sys.g(ieSpin,1:dgTensor) = g;
    end
  end
end

% Validate spin system
[System,err] = validatespinsys(Sys);
error(err);


% Build or load initial state
if isfield(System,'initState') && ~isempty(System.initState)
  % if some initial state was provided, this checks if the dimensions are
  % correct
  [a, b] = size(System.initState);
  if ischar(System.initState)
    error('String input for initial state not yet supported.')
  elseif System.nStates ~= a || System.nStates ~= b
    error('Initial state has to be a density matrix.')
  end
  Sigma = System.initState;
else
  % builds initial state, all electrons are -Sz, nuclei are not defined
  Sigma = -sop(System.Spins,'z1');
  for iElectron = 2 : System.nElectrons
    Sigma = Sigma - sop(System.Spins,['z' num2str(iElectron)]);
  end
end

% Build or load equilibtrium state - required for relaxation
if isfield(System,'eqState') && ~isempty(System.eqState)
  % if eqilibrium state was provided, this checks if the dimensions are
  % correct
  [a, b] = size(System.eqState);
  if ischar(System.eqState)
    error('String input for equilibrium state not yet supported.')
  elseif System.nStates ~= a || System.nStates ~= b
    error('Equilibrium state has to be a density matrix.')
  end
  Relaxation.equilibriumState = System.eqState;
else
  % initial state is copied from initial state
Relaxation.equilibriumState  = Sigma;
end


Ham = sham(System,B);

if isfield(Opt,'Relaxation') && ~isempty(Opt.Relaxation) && any(Opt.Relaxation)
  if isfield(System,'T1') || isfield(System,'T2')
    if ~isfield(System,'T1')
      System.T1 = 0;
    elseif ~isfield(System,'T2')
      System.T2 = 0;
    end
    Relaxation.Gamma = relaxationsuperoperator(System);
  else
    error('Relaxation was requested, but not relaxation times were provided.')
  end
end

if ~isfield(Det,'DetectionOperators') || isempty(Det.DetectionOperators)
  error('No detection operator provided.')
else
  Opt.DetectionOperators = Det.DetectionOperators;
end

nDetectionOperators = length(Opt.DetectionOperators);

DetectionOperators = cell(1,nDetectionOperators);

if ~isfield(Det,'FreqTranslation') || isempty(Det.FreqTranslation)
  Opt.FreqTranslation = [];
else
  Opt.FreqTranslation = zeros(1,nDetectionOperators);
  Opt.FreqTranslation(1:length(Det.FreqTranslation)) = Det.FreqTranslation;
end

Opt.DetectedEvents = zeros(1,length(Exp.t));

if ~isfield(Det,'Events') || isempty(Det.Events)
  Opt.DetectedEvents(end) = 1;
else
  Opt.DetectedEvents(1:length(Det.Events)) = Det.Events;
end

for iDetectionOperator = 1 : nDetectionOperators
  if ischar(Opt.DetectionOperators{iDetectionOperator})
    DetectionOperators{iDetectionOperator} = sop(System.Spins,Opt.DetectionOperators{iDetectionOperator});
  else
    DetectionOperators{iDetectionOperator} = Opt.DetectionOperators{iDetectionOperator};
  end
end


[Events, Vary] = sequencer(Exp,Opt);


% -------------------------------------------------------------------------
% Build the excitation operator for each pulse - if a custom excitation
% operator is provided in string form, sop will be called. If a matrix is
% provided, this matrix is used as excitation operator. If no custom
% excitation operators are given and no complex excitation was requested, 
% Sx is being used, for complex excitation Sx + Sy
% -------------------------------------------------------------------------

% initialize pulse counting
iPulse = 1;

% Loop over events and check if they are pulses
for iEvent = 1: length(Events)
  if strcmp(Events{iEvent}.type,'pulse')
    % Checks if user defined excitation operators were provided....
    if isfield(Opt,'ExcitationOperators') && ~isempty(Opt.ExcitationOperators) && (iPulse <= length(Opt.ExcitationOperators))
      % ... if yes, they are translated/verified and stored into the
      % current event ...
      if ischar(Opt.ExcitationOperators{iPulse})
        Events{iEvent}.xOp = sop(System.Spins,Opt.ExcitationOperators{iPulse});
      elseif any(size(Opt.ExcitationOperators{iPulse}) ~= size(Sigma))
        message = ['The excitation operator that you provided for pulse no. ' num2str(iPulse) ' does not have the same size as the density matrix.'];
        error(message)
      else
         Events{iEvent}.xOp = Opt.ExcitationOperators{iPulse};
      end
    else
      % ... if not, the Sx operator is formed for every spin in the
      % system...
      xOp = ['x' num2str(1)];
      Events{iEvent}.xOp = sop(System.Spins,xOp);
      for iSpin = 2 : length(System.Spins)
        xOp = ['x' num2str(iSpin)];
        Events{iEvent}.xOp = Events{iEvent}.xOp + sop(System.Spins,xOp);
      end
      % ... and if complex excitation is required, the Sy operator is added
      % on top
      if Events{iEvent}.ComplexExcitation
        for iSpin = 1 : length(System.Spins)
          yOp = ['y' num2str(iSpin)];
          Events{iEvent}.xOp= Events{iEvent}.xOp + sop(System.Spins,yOp);
        end
      end
    end
    % Increment pulse counter
    iPulse = iPulse + 1;
  end
end

% Calls the actual propagation engine
[TimeAxis, RawSignal, FinalState, StateTrajectories, NewEvents] = thyme(Sigma, Ham, DetectionOperators, Events, Relaxation, Vary);

% Signal postprocessing, such as down conversion and filtering
Signal = signalprocessing(TimeAxis,RawSignal,Opt);

end


