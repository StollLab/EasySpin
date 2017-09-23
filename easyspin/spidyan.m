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


% check for magnetic field
if length(Exp.Field) == 1
  B = [0 0 Exp.Field];
else
  B = Exp.Field;
end

if isfield(Opt,'FrameShift') && ~isempty(Opt.FrameShift)
  if isfield(Sys,'ZeemanFreq')
    Sys.ZeemanFreq =  Sys.ZeemanFreq - Opt.FrameShift;
  end
  
  if isfield(Sys,'g')  
    %%%%% Transfer from GHz (Opt.FrameShift) to MHz (Opt.FrameShift*1000)
    Sys.g = Sys.g - Opt.FrameShift*1000*1e9*planck/bmagn/Exp.Field(end);
  end
end

  %%%%% Transfer from GHz to MHz
Exp.Frequency = Exp.Frequency*1000;

% Translate Frequency to g values
if isfield(Sys,'ZeemanFreq')
  %%%%% Transfer from GHz to MHz
  Sys.ZeemanFreq = Sys.ZeemanFreq*1000;
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

if isfield(Sys,'ZeemanFreq')
  Sys = rmfield(Sys,'ZeemanFreq');
end

% This is spidyan specific
if ~isfield(Exp,'DetEvents') || isempty(Exp.DetEvents)
  Exp.DetEvents(length(Exp.t)) = 1;
else
  Exp.DetEvents(1:length(Exp.DetEvents)) = Exp.DetEvents;
end

[Events, Vary] = sequencer(Exp,Opt);


% Validate spin system
% use propagation setup instead here
[Sys,Sigma,DetOps,Events,Relaxation] = s_propagationsetup(Sys,Events,Opt);

Ham = sham(Sys,B);

% Spidyan specific
nDetOps = numel(DetOps);
if ~isfield(Opt,'FreqTranslation') || isempty(Opt.FreqTranslation)
  FreqTranslation = [];
else
  if isfield(Opt,'FrameShift') && ~isempty(Opt.FrameShift)
    Opt.FreqTranslation(Opt.FreqTranslation > 0) = Opt.FreqTranslation(Opt.FreqTranslation > 0) - Opt.FrameShift;
    Opt.FreqTranslation(Opt.FreqTranslation < 0) = Opt.FreqTranslation(Opt.FreqTranslation < 0) + Opt.FrameShift;
  end
  FreqTranslation = zeros(1,nDetOps);
  FreqTranslation(1:length(Opt.FreqTranslation)) = Opt.FreqTranslation;
end


% Calls the actual propagation engine
[TimeAxis, RawSignal, FinalState, StateTrajectories, NewEvents] = thyme(Sigma, Ham, DetOps, Events, Relaxation, Vary);

% Signal postprocessing, such as down conversion and filtering
Signal = signalprocessing(TimeAxis,RawSignal,DetOps,FreqTranslation);

% Signal = squeeze(Signal);

end


