% stochtraj  Generate stochastic rotational trajectories
%
%  [t,RTraj] = stochtraj(Par)
%  [t,RTraj,qTraj] = stochtraj(...)
%
%  Sys: stucture with system's dynamical parameters
%     .tcorr          correlation time (in seconds, 1 or 3 elements)
%     .lambda         ordering potential coefficient
%
%  Par: structure with simulation parameters
%     .dt             time step (in seconds)
%     .nSteps         number of time steps per simulation
%     .nTraj          number of trajectories
%     .alpha            Euler angle phi for starting orientation(s)
%     .beta          Euler angle theta for starting orientation(s)
%     .gamma            Euler angle chi for starting orientation(s)
%     .seed           seed the random number generator for reproducibility
%     .chkcon         if equal to 1, check for convergence in the angular
%                     distribution after every nSteps, and extend 
%                     propagation by an additional nSteps until it has
%                     converged
%     .Verbosity      0: no display, 1: show info
%
%  Output:
%     t              time points of the trajectory (in seconds)
%     RTraj              array of rotation matrices
%     qTraj          array of normalized quaternions

% Implementation based on 
%   Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%     http://dx.doi.org/10.1063/1.2908075

function varargout = stochtraj(Par)

if (nargin == 0), help(mfilename); return; end

if (nargin > 1), error('Too many input arguments.'); end

switch nargout
  case 0 % plotting
  case 2 % t,RTraj
  case 3 % t,RTraj,qTraj
  otherwise
    error('Incorrect number of output arguments.');
end

if ~isfield(Par,'Verbosity')
  Par.Verbosity = 0; % Log level
end

global EasySpinLogLevel;
EasySpinLogLevel = Par.Verbosity;

%========================================================================
% Dynamics and ordering potential
%========================================================================
% if no ordering potential coefficient is given, set it to 0
if ~isfield(Par,'lambda'), Par.lambda = 0; end

if numel(Par.lambda) > 1
  error('Only one orienting potential coefficient, c20, is currently implemented.')
end
Sim.lambda = Par.lambda;

% parse the dynamics parameter input using private function
if isfield(Par,'tcorr'), Dynamics.tcorr = Par.tcorr;
elseif isfield(Par,'Diff'), Dynamics.Diff = Par.Diff;
elseif isfield(Par,'logtcorr'), Dynamics.logtcorr = Par.logtcorr;
elseif isfield(Par,'logDiff'), Dynamics.logDiff = Par.logDiff;
else error('A rotational correlation time or diffusion rate is required.'); end

[Dynamics, err] = processdynamics(Dynamics);
error(err);
Sim.Diff = Dynamics.Diff';

%========================================================================
% Discrete Monte carlo settings
%========================================================================
if isfield(Par,'t')
  % time axis is given explicitly
  t = Par.t;
  Sim.nSteps = numel(t);
  Sim.dt = t(2) - t(1);
  if (abs(Sim.dt - max(t)/Sim.nSteps) > eps), error('t is not linearly spaced.'); end

elseif isfield(Par,'nSteps') && isfield(Par,'dt')
  % number of steps and time step are given
  Sim.dt = Par.dt;
  Sim.nSteps = Par.nSteps;
  t = linspace(0, Sim.nSteps*Sim.dt, Sim.nSteps);

elseif isfield(Par, 'nSteps') && isfield(Par, 'tmax')
  % number of steps and max time are given
  tmax = Par.tmax;
  Sim.nSteps = Par.nSteps;
  Sim.dt = tmax/Sim.nSteps;
  t = linspace(0, Sim.nSteps*Sim.dt, Sim.nSteps);

else
  error(['You must specify a time array, or a number of steps and ' ...
        'either a time step or tmax.'])
end

if isfield(Par,'seed')
  rng(Par.seed);
end

%========================================================================
% Grid and trajectory settings
%========================================================================
% If number of trajectories is not given, set it to 1
if ~isfield(Par, 'nTraj'), Par.nTraj = 1; end
Sim.nTraj = Par.nTraj;

% Get user-supplied starting angles
alpha = [];
beta = [];
gamma = [];
if isfield(Par,'alpha'), alpha = Par.alpha; end
if isfield(Par,'beta'), beta = Par.beta; end
if isfield(Par,'gamma'), gamma = Par.gamma; end

% Supplement starting angles if necessary
if isempty(alpha)
  alpha = rand(1,Sim.nTraj)*2*pi;
end
if isempty(beta)
  z = 2*rand(1,Sim.nTraj)-1;
  beta = acos(z);
end
if isempty(gamma)
  % Orienting potentials with m'=0 are independent of gamma, so we can set chi0=0
  gamma = zeros(1,Sim.nTraj);
end

% If only one starting angle and multiple trajectories, repeat the angle
if numel(beta) == 1 && Sim.nTraj > 1
  beta = repmat(beta,1,Sim.nTraj);
elseif numel(beta) ~= Sim.nTraj
  error('The number of starting angles must be equal to the number of trajectories.');
end
if numel(alpha) == 1 && Sim.nTraj > 1
  alpha = repmat(alpha,1,Sim.nTraj);
elseif numel(alpha) ~= Sim.nTraj
  error('The number of starting angles must be equal to the number of trajectories.');
end
if numel(gamma) == 1 && Sim.nTraj > 1
  gamma = repmat(gamma,1,Sim.nTraj);
elseif numel(gamma) ~= Sim.nTraj
  error('The number of starting angles must be equal to the number of trajectories.');
end

assert(numel(beta) == numel(alpha), 'Beta and alpha must be the same size.')
assert(numel(beta) == numel(gamma), 'Beta and gamma must be the same size.')

q0 = euler2quat(alpha,beta,gamma);

% Convert initial quaternion to a unitary 2x2 matrix for easier propagation
Q = zeros(2,2,Sim.nSteps,Sim.nTraj);
Q(1,1,1,:) =  q0(1,:) - 1i*q0(4,:);
Q(1,2,1,:) = -q0(3,:) - 1i*q0(2,:);
Q(2,1,1,:) =  q0(3,:) - 1i*q0(2,:);
Q(2,2,1,:) =  q0(1,:) + 1i*q0(4,:);

% if isfield(Par,'chkcon')
%   chkcon = Par.chkcon;
% else
%   chkcon = 0;
% end

% if isfield(Par,'tol')
%   if ~isfield(Par,'chkcon') || chkcon==0
%     error('A convergence tolerance was provided, but ''Par.chkcon'' was not equal to 1.')
%   end
%   tol = Par.tol;
% else
%   tol = 1e-3;
% end
  

%========================================================================
% Begin simulation
%========================================================================

converged = 0;
iter = 0;

logmsg(1,'-- Calculating stochastic trajectories -----------------------');

while ~converged

  % Pre-calculate angular steps due to random torques
  %   (Eq. 61 from reference, without factor of 1/2)
  Sim.randAngStep = bsxfun(@times, randn(3,Sim.nSteps-1,Sim.nTraj),...
                                   sqrt(2*Sim.Diff*Sim.dt));
    
  if ~isfield(Sim,'lambda')||Sim.lambda==0
    % No ordering potential -> only random torque present
    %  (Eqs. 48 and 61 in reference)
    Q = free_diff(Q, Sim);

  else
    % Ordering potential present -> include systematic torque
    %  (Eqs. 48 and 61 in reference)
    Q = aniso_diff(Q, Sim);
  
  end
  
  if iter==0
    % Convert 2x2 matrices to 4x1 quaternions
    qTraj = zeros(4,Sim.nSteps,Sim.nTraj);
    qTraj(1,:,:) = squeeze(real(Q(1,1,:,:)));
    qTraj(2,:,:) = squeeze(-imag(Q(1,2,:,:)));
    qTraj(3,:,:) = squeeze(real(Q(2,1,:,:)));
    qTraj(4,:,:) = squeeze(imag(Q(2,2,:,:)));

    % Convert to rotation matrices
    RTraj = quat2rotmat(qTraj);
  else
    % If not converged, then extend simulation along time axis
    qTemp(1,:,:) = squeeze(real(Q(1,1,:,:)));
    qTemp(2,:,:) = squeeze(-imag(Q(1,2,:,:)));
    qTemp(3,:,:) = squeeze(real(Q(2,1,:,:)));
    qTemp(4,:,:) = squeeze(imag(Q(2,2,:,:)));
    
    qTraj = cat(3,qTraj,qTemp);

    RTraj = cat(4,RTraj,quat2rotmat(qTemp));
  end

chkcon = 0;  

  if chkcon
    sprintf('Convergence tests not implemented yet!')
    converged = 1;
    [converged, ChiSquare] = check_convergence(RTraj, Sim.lambda, tol);
  else
    converged = 1;
  end
  
%   iter = iter + 1;
%   
%   if iter>10
%     logmsg(1,'Warning: convergence is very slow. Consider increasing\nlength of trajectories or increasing the tolerance value.')
%   end

end

% RTraj = quat2rotmat(qTraj);

logmsg(1,'Propagation was extended %d times before converging.', iter-1)

%==============================================================
%  Final processing
%==============================================================

switch nargout
  case 0 % Plot results
    maxTraj = 3;
    if Sim.nTraj>maxTraj
      error('Cannot plot more than %d trajectory.',maxTraj);
    end
    clf
    hold on
    for iTraj = 1:min(maxTraj,Sim.nTraj)
      x = squeeze(RTraj(1,3,:,iTraj));
      y = squeeze(RTraj(2,3,:,iTraj));
      z = squeeze(RTraj(3,3,:,iTraj));
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
  
  case 2  % Output rotation matrices only
    varargout = {t, RTraj};
    
  case 3  % Output rotation matrices and quaternions
    varargout = {t, RTraj, qTraj};
    
end

logmsg(1,'--------------------------------------------------------------');

clear global EasySpinLogLevel

return

%========================================================================
%========================================================================
function [Dyn,err] = processdynamics(D)%,FieldSweep)

Dyn = D;
err = '';

% diffusion tensor, correlation time
%------------------------------------------------------------------------
% convert everything (tcorr, logcorr, logDiff) to Diff
if isfield(Dyn,'Diff')
  % Diff given
elseif isfield(Dyn,'logDiff')
  Dyn.Diff = 10.^Dyn.logDiff;
elseif isfield(Dyn,'tcorr')
  Dyn.Diff = 1/6./Dyn.tcorr;
elseif isfield(Dyn,'logtcorr')
  if Dyn.logtcorr>=0, error('Sys.logtcorr must be negative.'); end
  Dyn.Diff = 1/6./10.^Dyn.logtcorr;
else
  err = sprintf('You must specify a rotational correlation time or a diffusion tensor\n(Sys.tcorr, Sys.logtcorr, Sys.Diff or Sys.logDiff).');
  return
end

if any(Dyn.Diff<0)
  error('Negative diffusion rate or correlation times are not possible.');
elseif any(Dyn.Diff>1e12)
  fprintf('Diffusion rate very fast. Simulation might not converge.\n');
elseif any(Dyn.Diff<1e3)
  fprintf('Diffusion rate very slow. Simulation might not converge.\n');
end

% expand to rhombic tensor
switch numel(Dyn.Diff)
  case 1, Dyn.Diff = Dyn.Diff([1 1 1]);
  case 2, Dyn.Diff = Dyn.Diff([1 1 2]);
  case 3, % Diff already rhombic
  otherwise
    err = 'Sys.Diff must have 1, 2 or 3 elements (isotropic, axial, rhombic).';
    return
end

end

%========================================================================

function Q = free_diff(Q, Sim)
% Propagate Q matrix assuming free, Brownian diffusion
randAngStep = Sim.randAngStep;
nSteps = Sim.nSteps;
nTraj = Sim.nTraj;

% Calculate size and normalized axis of angular step
theta = sqrt(sum(randAngStep.^2,1));
ux = randAngStep(1,:,:)./theta;
uy = randAngStep(2,:,:)./theta;
uz = randAngStep(3,:,:)./theta;

st = sin(theta/2);
ct = cos(theta/2);

% Pre-allocate array of propagators
% All propagators can be pre-calculated for Brownian diffusion
U = zeros(2,2,nSteps-1,nTraj);

U(1,1,:,:) = ct - 1i*uz.*st;
U(1,2,:,:) = -(uy + 1i*ux).*st;
U(2,1,:,:) = (uy - 1i*ux).*st;
U(2,2,:,:) = ct + 1i*uz.*st;


% Perform propagation
for iStep = 2:nSteps
  Q(:,:,iStep,:) = matmult(Q(:,:,iStep-1,:),U(:,:,iStep-1,:));
end

end

%========================================================================

function Q = aniso_diff(Q, Sim)
% Propagate Q matrix assuming Brownian diffusion in an orienting potential
randAngStep = Sim.randAngStep;
nSteps = Sim.nSteps;
nTraj = Sim.nTraj;
dt = Sim.dt;
Diff = Sim.Diff;
lambda = Sim.lambda;

for iStep=2:nSteps
  
  torque = anistorque(Q(:,:,iStep-1,:), lambda);
  if iter==0
    AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,iStep-1,:);
  else
    AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,iStep,:);
  end

  % Calculate size and normalized axis of angular step
  theta = sqrt(sum(AngStep.^2, 1));
  ux = AngStep(1,:,:)./theta;
  uy = AngStep(2,:,:)./theta;
  uz = AngStep(3,:,:)./theta;
  
  st = sin(theta/2);
  ct = cos(theta/2);

  % For anisotropic potentials, propagators must be calculated for each
  % timestep
  % Need '1' in the third dimension so that size(U) matches size(Q(:,:,iStep-1,:))
  U(1,1,1,:) = ct - 1i*uz.*st;
  U(1,2,1,:) = -(uy + 1i*ux).*st;
  U(2,1,1,:) = (uy - 1i*ux).*st;
  U(2,2,1,:) = ct + 1i*uz.*st;

  % Perform propagation
  Q(:,:,iStep,:) = matmult(Q(:,:,iStep-1,:), U);
end

end

%========================================================================

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

%========================================================================

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

%========================================================================
    
end