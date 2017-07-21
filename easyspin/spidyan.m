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

function [t, signal,state,sigmas,Eventsnew] = spidyan(Sys,Exp,Opt)

% if (nargin==0), help(mfilename); return; end

% Get time for performance report at the end.
StartTime = clock;

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


% Build the Spin System Here
options.silent = 1;
options.relaxation = 1;

sqn = 1/2;
% Det = {spops(sqn,'x') spops(sqn,'z') spops(sqn,'p')/2};
Det = {spops(1/2,'z')};
% 
system.sqn = sqn;
system.interactions = {1,0,'z','e',1.5};
system.T1 = 1;
system.T2 = 5;
[system,Sigma] = setup(system,options);
Ham = system.ham*1000/2/pi;
Relaxation.equilibriumState = system.eq;
Relaxation.Gamma = system.gamma;

[Events, Vary] = sequencer(Exp,Opt);

% %And the vary structure
% Vary.Events = {[1] 2};
% Vary.IQs{1} = {Events{1}.IQ Events{1}.IQ/2};
% % Vary.IQs{3} = {IQ/4 IQ};
% Vary.ts{1} = {Events{1}.t Events{1}.t};
% % Vary.ts{3} = {t t};
% Vary.Points = [2 2];
% 
% Vary.ts{2} = {0:Exp.TimeStep:0.2; 0:Exp.TimeStep:0.4};

[t, signal,state,sigmas,Eventsnew]=evolve2(Sigma, Ham, Det, Events, Relaxation, Vary);