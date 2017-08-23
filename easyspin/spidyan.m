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

function [TimeAxis, Signal,FinalState,StateTrajectories,NewEvents] = spidyan(Sys,Exp,Opt)

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



if length(Exp.B) == 1
  B = [0 0 Exp.B];
else
  B = Exp.B;
end

if isfield(Sys,'ResonanceFrequency')
  %Translate Frequency to g values
  if isfield(Sys,'g')
    [~, dgTensor] = size(Sys.g);
  else
    dgTensor = 1;
  end
  ReferenceSpin = struct('S',1/2,'g',2);
  RefZeeman = zeeman(ReferenceSpin,B);
  RefEigenValues = eig(RefZeeman);
  RefOmega = RefEigenValues(2)-RefEigenValues(1);
  
  for ieSpin = 1 : length(Sys.ResonanceFrequency)
    if Sys.ResonanceFrequency(ieSpin) ~= 0
      Adaptedg = Sys.ResonanceFrequency(ieSpin)*2/RefOmega;
      Sys.g(ieSpin,1:dgTensor) = Adaptedg;
    end
  end
end

% Validate spin system
[System,err] = validatespinsys(Sys);
error(err);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIDYAN commands
options.silent = 1;
options.relaxation = 1;
sqn = System.Spins;
system.sqn = sqn;
system.interactions = {1,0,'z','e',1.5};
system.T1 = 1;
system.T2 = 5;
[system,Sigma] = setup(system,options);
% Ham = system.ham*1000/2/pi;
Relaxation.equilibriumState = system.eq;
Gamma = system.gamma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nDetectionOperators = length(Opt.DetectionOperators);

DetectionOperators = cell(1,nDetectionOperators);

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
[TimeAxis, RawSignal, FinalState, StateTrajectories, NewEvents]=thyme(Sigma, Ham, DetectionOperators, Events, Relaxation, Vary);

% Signal postprocessing, such as down conversion and filtering
Signal = signalprocessing(TimeAxis,RawSignal,Opt);

end


