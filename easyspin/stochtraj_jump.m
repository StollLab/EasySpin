% stochtraj_jump  Generate stochastic trajectories of Markovian jumps using
%                 kinetic Monte Carlo
%
%   [t,RTraj] = stochtraj_jump(Sys)
%   ... = stochtraj_jump(Sys,Par)
%   ... = stochtraj_jump(Sys,Par,Opt)
%   [t,RTraj,qTraj] = stochtraj_jump(...)
%   [t,RTraj,qTraj,stateTraj] = stochtraj_jump(...)
%   [t,stateTraj] = stochtraj_jump(...)
%
%   Sys: stucture with system's dynamical parameters
%
%     TransRates     numeric, size = (nStates,nStates)
%                    transition rate matrix describing inter-state dynamics
%                    for kinetic Monte Carlo simulations
%
%     TransProb      numeric, size = (nStates,nStates)
%                    transition probability matrix describing inter-state 
%                    dynamics for kinetic Monte Carlo simulations, note
%                    that a time step must be given to use Sys.TransProb
%                    (alternative input to TransRates; ignored if TransRates
%                    is given)
%
%     Orientations   numeric, size = (3,nStates)
%                    Euler angles for each state's orientation
%
%
%   Par: simulation parameters for Monte Carlo integrator
%
%     nTraj          integer
%                    number of trajectories
%
%     StatesStart    numeric, size = (1,1) or (nTraj,1)
%                    starting states of trajectories
%
%         When specifying the simulation time provide one of the following
%         combinations:
%         Precedence: nSteps,dt > tMax,nSteps
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
%     statesOnly        1 or 0
%                       specify whether or not to use a discrete model purely
%                       to generate states and not quaternion orientations
%
% %     checkConvergence  if equal to 1, after the first nSteps of the 
% %                       trajectories are calculated, both inter- and intra-
% %                       trajectory convergence is checked using the Gelman-
% %                       Rubin R statistic such that R<1.1, and if this 
% %                       condition is not satisfied, then propagation will be 
% %                       extended by either a length of time equal to the 
% %                       average of tcorr or by 20% more time steps, whichever 
% %                       is larger
%
%     Verbosity         0: no display, 1: show info
%
%
%   Output:
%
%     t              matrix, size = (nSteps,1) 
%                    time points of the trajectory (in seconds)
%
%     RTraj          4D array, size = (3,3,nSteps,nTraj)
%                    trajectories of rotation matrices
%
%     qTraj          3D array, size = (4,nSteps,nTraj)
%                    trajectories of normalized quaternions
%
%     stateTraj      2D array, size = (nSteps,nTraj)
%                    trajectories of states

function varargout = stochtraj_jump(Sys,Par,Opt)

% Preprocessing
% -------------------------------------------------------------------------

% Only Sys needs to be given to run stochtraj_jump, so if Par is 
% not given, initialize it here
switch nargin
  case 0
    help(mfilename); return;
  case 1
    % only Sys is given
    Par = struct;
    Opt = struct;
  case 2
    % Opt is not given
    Opt = struct;
  case 3
    % do nothing
  otherwise
    error('Too many input arguments.')
end

switch nargout
  case 0 % plotting
  case 2 % t,RTraj; t,stateTraj
  case 3 % t,RTraj,qTraj
  case 4 % t,RTraj,qTraj,stateTraj
  otherwise
    error('Incorrect number of output arguments.');
end

% Check Opt
% -------------------------------------------------------------------------

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

if ~isfield(Opt,'statesOnly')
  Opt.statesOnly = false || nargout==0;
end

global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;


% Check dynamics
% -------------------------------------------------------------------------

if isfield(Sys,'Potential')
  warning('Sys.Potential is given. This field is not used by stochtraj_jump.');
end

if ~isfield(Par,'dt')
  error('The time step Par.dt must be specified.')
end

if isfield(Sys,'TransRates')
  
  TRM = Sys.TransRates;
  
  % error checking
  if ~isnumeric(TRM) || ~ismatrix(TRM) || size(TRM,1)~=size(TRM,2)
    error('Sys.TransRates must be a square matrix.')
  end
  any2d = @(M)any(M(:));
  if any(diag(TRM)>0) || any2d(TRM-diag(diag(TRM))<0)
    error('Sys.TransRates must contain strictly negative diagonal and positive off-diagonal entries.')
  end
  if any(abs(sum(TRM,2)/max(abs(TRM(:))))>1e-13)
    error('In Sys.TransRates, the sum over each row must be zero.')
  end
  
  TPM = expm(Par.dt*TRM);

  % get the relaxation times
  D = eig(TRM);
  idx = abs(D)/max(abs(D))>1e-11;
  tcorr = -1./D(idx);
  
elseif isfield(Sys,'TransProb')
  
  TPM = Sys.TransProb;
  
  % error checking
  if ~isnumeric(TPM) || ~ismatrix(TPM) || size(TPM,1)~=size(TPM,2)
    error('Sys.TransProb must be a square matrix.')
  end
  if any(TPM(:)<0)
    error('All elements in Sys.TransProp must positive (or zero).');
  end
  if any(abs(nonzeros(sum(TPM,2))-1)>1e-12)
    error('The rows of Sys.TransProb must sum to 1.')
  end
  
else
  
  error(['A transition rate matrix (Sys.TransRates) or a transition probability matrix (Sys.TransProb) ',... 
         'is required.'])
       
end

nStates = size(TPM,1);


% set kinetic Monte Carlo cumulative transition probability matrix
% cumulTPM = cumsum(TPM,1);
cumulTPM = cumsum(TPM,2);
  
if ~Opt.statesOnly
  if isfield(Sys,'Orientations')
    Orientations = Sys.Orientations;
    if isequal(size(Orientations),[3,nStates])
      % do nothing
    elseif isequal(size(Orientations),[nStates,3])
      % transpose to size (3,nStates)
      Orientations = Orientations.';
    else
      error('The size of Sys.Orientations must be (3,nStates) or (nStates,3).')
    end
  else
    error(['Orientations for the %d states are required. Give a %dx%d array in Sys.Orientations.',...
      '\nIf only states are desired, use Opt.statesOnly = 1.'],nStates,nStates,3);
  end
end


% Discrete Monte Carlo settings (Par)
% -------------------------------------------------------------------------

% set integration time step and number of simulation steps
  
if isfield(Par,'nSteps') && isfield(Par,'dt')
  % number of steps and time step are given
  dt = Par.dt;
  nSteps = Par.nSteps;

elseif isfield(Par,'tMax') && isfield(Par,'dt')
  dt = Par.dt;
  if (dt<=0)
    error('Par.dt must be positive.');
  end
  nSteps = ceil(Par.tMax/dt);

else
  % only time step is given
  dt = Par.dt;
  if isfield(Sys,'TransRates')
    nSteps = ceil(200*max(tcorr)/dt);
  elseif isfield(Sys,'TransProb')
    nSteps = 500;
  end
  logmsg(0,'-- Number of time steps not given. Using %d steps.', nSteps);
  
end

Par.nSteps = nSteps;
Par.dt = dt;

% set integrator type
if ~isfield(Par,'Integrator')
  % default Monte Carlo integrator is Euler-Maruyama
  Par.Integrator = 'Euler-Maruyama';
end

if ~strcmp(Par.Integrator,'Euler-Maruyama')&&~strcmp(Par.Integrator,'Leimkuhler-Matthews')
  error('Input for integrator method not recognized.')
end


% Grid and trajectory settings
% -------------------------------------------------------------------------

% If number of trajectories is not given, set it to 1
if ~isfield(Par,'nTraj'), Par.nTraj = 1; end

% Get user-supplied starting states
if isfield(Par,'StatesStart')
  StatesStart = Par.StatesStart;
  if ~isvector(StatesStart) || numel(StatesStart)~=Par.nTraj
    error('Par.StatesStart should be of size (1,1) or (1,Par.nTraj).')
  end
  if any(StatesStart<1) || any(StatesStart>nStates)
    error(['Each entry in Par.StatesStart needs to be equal to an integer ',...
           'between 1 and nStates.'])
  end
  StatesStart = StatesStart(:);
else
  StatesStart = randi(nStates,1,Par.nTraj);
end

if isfield(Opt,'checkConvergence')
  checkConvergence = Opt.chkcon;
  if checkConvergence && Par.nTraj==1
    error(['Checking convergence of a single trajectory using the ',...
           'R statistic is not supported.\n'])
  end
else
  checkConvergence = false;
end

if ~Opt.statesOnly
  qStates = euler2quat(Orientations);
  % initialize quaternion trajectories and their starting orientations
  qTraj = zeros(4,Par.nSteps,Par.nTraj);
  for iTraj = 1:Par.nTraj
    qTraj(:,1,iTraj) = qStates(:,StatesStart(iTraj));
  end
end
 

% Simulation
% -------------------------------------------------------------------------
logmsg(2,'-- Calculating stochastic trajectories -----------------------');

nTraj = Par.nTraj;
nSteps = Par.nSteps;

stateTraj = zeros(nSteps,nTraj);
stateTraj(1,:) = StatesStart;

statesOnly = Opt.statesOnly;
for iTraj = 1:nTraj
  u = rand(1,nSteps);
  state = stateTraj(1,iTraj);
  for iStep = 2:nSteps
    state = find(cumulTPM(state,:)>u(iStep-1),1);
    stateTraj(iStep,iTraj) = state;
    if ~statesOnly
      qTraj(:,iStep,iTraj) = qStates(:,state);
    end
  end
end

t = linspace(0, nSteps*Par.dt, nSteps).';

logmsg(2,'-- Propagation finished --------------------------------------');
logmsg(2,'--------------------------------------------------------------');


% Final processing
% -------------------------------------------------------------------------

switch nargout
  case 0 % Plot results
    maxTraj = 3;
    
    clf
    
    nPlotTraj = min(nTraj,maxTraj);
    for iTraj = 1:nPlotTraj
      subplot(nPlotTraj,1,iTraj)
      hold on
      s = stateTraj(:,iTraj);
      
      hl = plot(t,s);
      set(hl,'Color',[1 1 1]*0.3);
      for iState = 1:nStates
        idx = s==iState;
        h = plot(t(idx),iState*ones(1,sum(idx)),'o');
      end
    
      xlabel('time (s)');
      ylabel('state');
      set(gca,'YTick',1:nStates);
      axis tight
      stateLims = [1 nStates] + [-1 1]*0.5;
      ylim(stateLims);
      box on
      title(sprintf('Trajectory %d (%d states, %d steps)',iTraj,nStates,length(s)));
      
      histdata = zeros(1,nStates);
      for iState = 1:nStates
        histdata(iState) = sum(s==iState);
      end
      
      p0 = get(gca,'Position');
      p0(3) = 0.7;
      set(gca,'Position',p0);
      p1 = p0;
      p1(1) = p0(1)+p0(3);
      p1(3) = (1-p0(3)-p0(1))-0.07;
      ax = axes('Position',p1);
      barh(1:nStates,histdata,0.3);
      ylim(stateLims)
      set(ax,'XAxisLocation','top');
      set(ax,'YAxisLocation','right');
      ylabel(ax,'histogram');
      xlabel(ax,'count')
      set(ax,'YTick',[]);
      %ax.YLabel.Rotation = -90;
    end

  case 2  % Output rotation matrix or state trajectories
    if Opt.statesOnly
      varargout = {t, stateTraj};
    else
      RTraj = quat2rotmat(qTraj);
      varargout = {t, RTraj};
    end
  
  case 3  % Output both rotation matrix and quaternion trajectories
    RTraj = quat2rotmat(qTraj);
    varargout = {t, RTraj, qTraj};
    
  case 4  % Output rotation matrix, quaternion and state trajectories
    RTraj = quat2rotmat(qTraj);
    varargout = {t, RTraj, qTraj, stateTraj};

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