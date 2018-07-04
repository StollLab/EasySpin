% stochtraj_jump  Generate stochastic trajectories of Markovian jumps using
%                 kinetic Monte Caro
%
%   [t,qTraj] = stochtraj_jump(Sys)
%   [t,qTraj] = stochtraj_jump(Sys,Par)
%   [t,qTraj] = stochtraj_jump(Sys,Par,Opt)
%   [t,stateTraj] = stochtraj_jump(...)
%   [t,qTraj,stateTraj] = stochtraj_jump(...)
%
%   Sys: stucture with system's dynamical parameters
%
%     TransRates     numeric, size = (nStates,nStates)
%                    transition rate matrix describing inter-state dynamics
%                    for kinetic Monte Carlo simulations
%
%     TransProb      numeric, size = (nStates,nStates)
%                    transition probability matrix describing inter-state 
%                    dynamics for kinetic Monte Carlo simulations
%
%     States         numeric, size = (3,nStates)
%                    Euler angles for each state's orientation
%
%
%
%   Par: simulation parameters for Monte Carlo integrator
%
%     nTraj          integer
%                    number of trajectories
%
%     States0        numeric, size = (1,1) or (nTraj,1)
%                    starting states of Markov chains
%
%         When specifying the simulation time provide one of the following
%         combinations:
%         Precedence: t > nSteps,dt > tMax,nSteps
%             
%     t              numeric
%                    array of time points  (in seconds)
%
%     tMax           double
%                    total time of simulation (in seconds)
%
%     dt             double
%                    time step (in seconds)
%
%     nSteps         integer
%                    number of time steps in simulation
%
%
%   Opt: simulation options
%
%     statesOnly     1 or 0
%                    specify whether or not to use a discrete model purely
%                    to generate states and not quaternion orientations
%
%     chkcon         if equal to 1, after the first nSteps of the 
%                    trajectories are calculated, both inter- and intra-
%                    trajectory convergence is checked using the Gelman-
%                    Rubin R statistic such that R<1.1, and if this 
%                    condition is not satisfied, then propagation will be 
%                    extended by either a length of time equal to the 
%                    average of tcorr or by 20% more time steps, whichever 
%                    is larger
%
%     Verbosity      0: no display, 1: show info
%
%
%   Output:
%
%     t              matrix, size = (nSteps,1) 
%                    time points of the trajectory (in seconds)
%
%     qTraj          3D array, size = (4,nTraj,nSteps)
%                    trajectories of normalized quaternions
%
%     stateTraj      3D array, size = (nStates,nTraj,nSteps)
%                    trajectories of states

function varargout = stochtraj_jump(Sys,Par,Opt)

% Preprocessing
% -------------------------------------------------------------------------

% Only Sys needs to be given to run stochtraj_jump properly, so if Par is 
% not given, initialize it here
switch nargin
  case 0
    help(mfilename); return;
  case 1
    % only Sys is given
    Par = struct('unused',NaN);
    Opt = struct('unused',NaN);
%    skipchk = 0;
  case 2
    % Opt is not given
    Opt = struct('unused',NaN);
%    skipchk = 0;
  case 3
    % do nothing
%    skipchk = 0;
%   case 4
%     % 4th argument is used as a flag to skip error checking
%     if ischar(flag) && strcmp(flag,'skipchk')
%       skipchk = 1;
%     else
%       error('Fourth argument is for internal use only. Please remove it.'); 
%     end
  otherwise
    error('Too many input arguments.')
end

switch nargout
  case 0 % plotting
  case 2 % t,qTraj
%   case 2 % t,qTraj or stateTraj
  case 3 % t,qTraj,stateTraj
  otherwise
    error('Incorrect number of output arguments.');
end

% Check Opt
% -------------------------------------------------------------------------

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

if ~isfield(Opt,'statesOnly')
  Opt.statesOnly = 0;
end


global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;


% Check dynamics
% -------------------------------------------------------------------------

if isfield(Sys,'TransRates')
  TRM = Sys.TransRates;
  if ~isnumeric(TRM) || ~ismatrix(TRM) || (size(TRM,1)~=size(TRM,2))
    error('TransRates must be a square matrix.')
  end
  nStates = size(TRM,1);
  diagRates = diag(diag(TRM));
  if any(diag(TRM)>0) || any((TRM(:)-diagRates(:))<0)
    error('TransRates must contain strictly negative diagonal and positive off-diagonal entries.')
  end
  if any(abs(sum(TRM,2)/norm(TRM))>1e-13)
    error("In the TransRates matrix, the sum of each row's off-diagonal elements must equal the negative of the diagonal element.")
  end
elseif isfield(Sys,'TransProb')
  TPM = Sys.TransProb;
  if ~isnumeric(TPM) || ~ismatrix(TPM) || (size(TPM,1)~=size(TPM,2))
    error('TransProb must be a square matrix.')
  end
  nStates = size(TPM,1);
else
  error(['A transition rate matrix or a transition probability matrix ',... 
         'is required a jump simulation.'])
end
  
if ~Opt.statesOnly
  if isfield(Sys,'States')
    States = Sys.States;
    if size(States,1)~=3 || size(States,2)~=nStates
      error(['The size of States must be (3,nStates), with the size of the ' ...
             'second dimension equal to the number of rows (and columns) of TransProb.'])
    end
  else
    error('A set of States is required for a jump simulation.')
  end
end


% Discrete Monte carlo settings (Par)
% -------------------------------------------------------------------------

% set integration time step and number of simulation steps
if isfield(Par,'t')
  % time axis is given explicitly
  t = Par.t;
  nSteps = numel(t);
  dt = t(2) - t(1);
  if (abs(dt - max(t)/nSteps) > eps), error('t is not linearly spaced.'); end
  
elseif isfield(Par,'nSteps') && isfield(Par,'dt')
  % number of steps and time step are given
  dt = Par.dt;
  nSteps = Par.nSteps;

elseif isfield(Par,'tMax') && isfield(Par,'nSteps')
  % number of steps and max time are given
  tMax = Par.tMax;
  nSteps = Par.nSteps;
  dt = tMax/Sim.nSteps;

else
%   error(['You must specify a time array, or a number of steps and ' ...
%         'either a time step or tmax.'])
  logmsg(1,'-- No time step given. Par.dt set to Par.tcorr/10: %0.5g s.', 1/6/mean(Sim.Diff));
  dt = 1/6/mean(Sim.Diff)/10;
  if ~isfield(Par, 'nSteps')
    logmsg(1,'-- Number of time steps not given. Par.nSteps set to 200e-9/Par.dt: %d.', ceil(200e-9/dt));
    nSteps = ceil(200e-9/dt);
  else
    nSteps = Par.nSteps;
  end
end

Sim.nSteps = nSteps;
Sim.dt = dt;

% set kinetic Monte Carlo cumulative transition probability matrix
if isfield(Sys,'TransRates')
  TPM = expm(Sim.dt*TRM);
end
cumulTPM = cumsum(TPM,1);
if any(abs(1-cumulTPM(end,:))>1e-13)
  error('The columns of cumulTPM = sum(TPM,1) must sum to 1.')
end

% set integrator type
if ~isfield(Par,'Integrator')
  % default Monte Carlo integrator is Euler-Maruyama
  Sim.Integrator = 'Euler-Maruyama';
else
  if ~strcmp(Par.Integrator,'Euler-Maruyama')&&~strcmp(Par.Integrator,'Leimkuhler-Matthews')
    error('Input for integrator method not recognized.')
  end
  Sim.Integrator = Par.Integrator;
end


% Grid and trajectory settings
% -------------------------------------------------------------------------

% If number of trajectories is not given, set it to 1
if ~isfield(Par, 'nTraj'), Par.nTraj = 1; end
Sim.nTraj = Par.nTraj;

% Get user-supplied starting states
if isfield(Par,'States0')
  States0 = Par.States0;
  if ~all(size(States0)==[1,1])&&~all(size(States0)==[1,Par.nTraj])
    error('States0 should be of size (1,1) or (1,Par.nTraj).')
  end
  if any(States0<1) || any(States0>nStates)
    error(['Each entry in States0 needs to be equal to an integer ',...
           'within the range [1,nStates].\n'])
  end
else
  States0 = randi(nStates,1,Par.nTraj);
end

if isfield(Opt,'chkcon')
  chkcon = Opt.chkcon;
  if chkcon==1 && Sim.nTraj==1
    error(['Checking convergence of a single trajectory using the ',...
           'R statistic is not supported.\n'])
  end
else
  chkcon = 0;
end

if ~Opt.statesOnly
  States = euler2quat(States);
  % initialize quaternion trajectories and their starting orientations
  qTraj = zeros(4,Sim.nTraj,Sim.nSteps);
  for iTraj = 1:Sim.nTraj
    qTraj(:,iTraj,1) = States(:,States0(iTraj));
  end
end

% if isfield(Par,'tol')
%   if ~isfield(Par,'chkcon') || chkcon==0
%     error('A convergence tolerance was provided, but ''Par.chkcon'' was not equal to 1.')
%   end
%   tol = Par.tol;
% else
%   tol = 1e-3;
% end
 

% Simulation
% -------------------------------------------------------------------------

converged = 0;

% if chkcon=1, then we need to keep track of how many iterations have been 
% completed (i.e. the number of instances of propagation in time by nSteps)
iter = 1;

logmsg(2,'-- Calculating stochastic trajectories -----------------------');

nTraj = Sim.nTraj;
nSteps = Sim.nSteps;

stateTraj = zeros(nTraj,nSteps);
stateTraj(:,1) = States0.';
u = rand(nTraj,nSteps);

for iTraj = 1:nTraj
  stateNew = stateTraj(iTraj,1);
  for iStep = 2:nSteps
    stateLast = stateNew;
    uLast = u(iTraj,iStep-1);
    stateNew = find(cumulTPM(:,stateLast)>uLast,1);
    stateTraj(iTraj,iStep) = stateNew;
    if ~Opt.statesOnly
      qTraj(:,iTraj,iStep) = States(:,stateNew);
    end
  end
end
totSteps = size(stateTraj,2);

t = linspace(0, totSteps*Sim.dt, totSteps).';

logmsg(2,'-- Propagation finished --------------------------------------');
logmsg(2,'--------------------------------------------------------------');


% Final processing
% -------------------------------------------------------------------------

switch nargout
  case 0 % Plot results
    maxTraj = 3;
    if Sim.nTraj>maxTraj
      error('Cannot plot more than %d trajectories.',maxTraj);
    end
    RTraj = quat2rotmat(q);
    clf
    hold on
    for iTraj = 1:min(maxTraj,Sim.nTraj)
      x = squeeze(RTraj(1,3,iTraj,:));
      y = squeeze(RTraj(2,3,iTraj,:));
      z = squeeze(RTraj(3,3,iTraj,:));
      plot3(x,y,z);
    end
    axis equal
    axlim = 1.2;
    xlim([-1 1]*axlim);
    ylim([-1 1]*axlim);
    zlim([-1 1]*axlim);
    ax = gca;
    ax.Box = 'on';
    ax.BoxStyle = 'full';
    view([-32 32]);
    xlabel('x');
    ylabel('y');
    zlabel('z');

  case 2  % Output quaternion trajectories
    if Opt.statesOnly
      varargout = {t, stateTraj};
    else
      varargout = {t, qTraj};
    end
  
  case 3  % Output both quaternion and state trajectories
    varargout = {t, qTraj, stateTraj};

end

clear global EasySpinLogLevel


end

% Helper functions
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------

% function varargout = acorr_convergence(RTraj, tol)
% % Calculate angular histogram of trajectories and compare with analytic
% % expression
% 
% nBins = 50;  % There might be a better value to choose here
%   
% VecTraj = squeeze(RTraj(:, 3, :, :));
% 
% totTraj = size(RTraj,4);
% 
% AutoCorrFFT = zeros(nSteps, nTraj);
% 
% for k = 1:totTraj
%   AutoCorrFFT(:, k) = autocorrfft(VecTraj(:, :, k).^2);
% end
% 
% AutoCorrFFT = sum(AutoCorrFFT, 2)'/totTraj;
% 
% converged = ChiSquare < tol;
% 
% varargout = {converged, ChiSquare};
% 
% end

% -------------------------------------------------------------------------

% function varargout = hist_convergence(RTraj, lambda, tol)
% % Calculate angular histogram of trajectories and compare with analytic
% % expression
% 
% nBins = 50;  % There might be a better value to choose here
%   
% VecTraj = squeeze(RTraj(:, 3, :, :));
% 
% totTraj = size(RTraj,4);
% 
% bins = linspace(0, pi, nBins)';
% ThetaHist = zeros(nBins, totTraj);
% 
% for k = 1:totTraj
%   ThetaHist(:, k) = hist(acos(VecTraj(3, :, k)), bins);
% end
% 
% ThetaHist = sum(ThetaHist, 2);
% ThetaHist = ThetaHist/sum(ThetaHist);
% 
% BoltzDist = exp(lambda*(1.5*cos(bins).^2 - 0.5));
% BoltzInt = sum(BoltzDist.*sin(bins));
% BoltzDist = BoltzDist.*sin(bins)./BoltzInt;
% 
% ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist);
% 
% converged = ChiSquare < tol;
% 
% varargout = {converged, ChiSquare};
% 
% end