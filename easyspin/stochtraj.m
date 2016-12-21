% stochtraj  Generate stochastic rotational trajectories
%
%  [t,RTraj] = stochtraj(Sys)
%  [t,RTraj] = stochtraj(Sys,Par)
%  [t,RTraj,qTraj] = stochtraj(...)
%
%  Sys: stucture with system's dynamical parameters
%     tcorr          correlation time (in seconds, 1 or 3 elements)
%
%     Coefs          array of coefficients
%
%     LMK            quantum numbers L, M, and K
%
%
%  Par: structure with simulation parameters
%     dt             time step (in seconds)
%
%     nSteps         number of time steps per simulation
%
%     nTraj          number of trajectories
%
%     alpha          Euler angle alpha for starting orientation(s)
%
%     beta           Euler angle beta for starting orientation(s)
%
%     gamma          Euler angle gamma for starting orientation(s)
%
%     seed           seed the random number generator for reproducibility
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
%  Output:
%     t              array, time points of the trajectory (in seconds)
%
%     RTraj          3x3x... array, rotation matrices
%
%     qTraj          4x... array, normalized quaternions

% Implementation based on 
%   Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%     http://dx.doi.org/10.1063/1.2908075

function varargout = stochtraj(Sys, Par)
%% Preprocessing

if (nargin == 0), help(mfilename); return; end

% Only Sys needs to be given to run stochtraj properly, so if Par is not 
% given, initialize it here
if (nargin == 1), Par.Verbosity = 0;
elseif (nargin > 2), error('Too many input arguments.'); end

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

%% Dynamics and ordering potential
%========================================================================

% if no ordering potential coefficient is given, set it to 0
if ~isfield(Sys,'Coefs') && ~isfield(Sys,'LMK')
  Sys.Coefs = [];
  Sys.LMK = [];
elseif xor(isfield(Sys,'Coefs'), isfield(Sys,'LMK'))
  error('Both ordering coefficients and LMK are required for an ordering potential.')
end

% if numel(Sys.lambda) > 1
%   error('Only one orienting potential coefficient, c20, is currently implemented.')
% end
Sim.Coefs = Sys.Coefs;
Sim.LMK = Sys.LMK;

% parse the dynamics parameter input using private function
if isfield(Sys,'tcorr'), Dynamics.tcorr = Sys.tcorr;
elseif isfield(Sys,'Diff'), Dynamics.Diff = Sys.Diff;
elseif isfield(Sys,'logtcorr'), Dynamics.logtcorr = Sys.logtcorr;
elseif isfield(Sys,'logDiff'), Dynamics.logDiff = Sys.logDiff;
else, error('A rotational correlation time or diffusion rate is required.'); end

[Dynamics, err] = processdynamics(Dynamics);
error(err);
Sim.Diff = Dynamics.Diff';
tcorrAvg = 1/6/mean(Sim.Diff);

%% Discrete Monte carlo settings
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

elseif isfield(Par, 'nSteps') && isfield(Par, 'tmax')
  % number of steps and max time are given
  tmax = Par.tmax;
  Sim.nSteps = Par.nSteps;
  Sim.dt = tmax/Sim.nSteps;

else
%   error(['You must specify a time array, or a number of steps and ' ...
%         'either a time step or tmax.'])
  logmsg(1,'-- No time step given. Par.dt set to Par.tcorr/10: %0.5g s.', 1/6/mean(Sim.Diff));
  Sim.dt = 1/6/mean(Sim.Diff)/10;
  if ~isfield(Par, 'nSteps')
    logmsg(1,'-- Number of time steps not given. Par.nSteps set to 200e-9/Par.dt: %d.', ceil(200e-9/Sim.dt));
    Sim.nSteps = ceil(200e-9/Sim.dt);
  else
    Sim.nSteps = Par.nSteps;
  end
end

if isfield(Par,'seed')
  rng(Par.seed);
end

%% Grid and trajectory settings
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

if isfield(Par,'chkcon')
  chkcon = Par.chkcon;
  if chkcon==1 && Sim.nTraj==1
    error('Checking convergence of a single trajectory using the R statistic is not supported.\n')
  end
else
  chkcon = 0;
end


q0 = euler2quat(alpha,beta,gamma);
qTraj = zeros(4,Sim.nTraj,Sim.nSteps);
qTraj(:,:,1) = q0;

% Convert initial quaternion to a unitary 2x2 matrix for easier propagation
% Q = zeros(2,2,Sim.nTraj,Sim.nSteps);
% Q(1,1,:,1) =  q0(1,:) - 1i*q0(4,:);
% Q(1,2,:,1) = -q0(3,:) - 1i*q0(2,:);
% Q(2,1,:,1) =  q0(3,:) - 1i*q0(2,:);
% Q(2,2,:,1) =  q0(1,:) + 1i*q0(4,:);

% if isfield(Par,'tol')
%   if ~isfield(Par,'chkcon') || chkcon==0
%     error('A convergence tolerance was provided, but ''Par.chkcon'' was not equal to 1.')
%   end
%   tol = Par.tol;
% else
%   tol = 1e-3;
% end
  
%% Simulation
%========================================================================

converged = 0;
iter = 0;

logmsg(1,'-- Calculating stochastic trajectories -----------------------');

while ~converged

  if iter==0
    %  Pre-calculate angular steps due to random torques
    %  (Eq. 61 from reference, without factor of 1/2)
    Sim.randAngStep = bsxfun(@times, randn(3,Sim.nTraj,Sim.nSteps),...
                                   sqrt(2*Sim.Diff*Sim.dt));

    %  Perform stochastic simulations
    %  (Eqs. 48 and 61 in reference)
%     Q = propagate(Q, Sim, iter);
    qTraj = propagate(qTraj, Sim, iter);
    
%     % Convert 2x2 matrices to 4x1 quaternions
%     qTraj = zeros(4,Sim.nTraj,Sim.nSteps);
%     qTraj(1,:,:) = squeeze(real(Q(1,1,:,:)));
%     qTraj(2,:,:) = squeeze(-imag(Q(1,2,:,:)));
%     qTraj(3,:,:) = squeeze(real(Q(2,1,:,:)));
%     qTraj(4,:,:) = squeeze(imag(Q(2,2,:,:)));

  else
    logmsg(1,'-- Convergence not obtained -------------------------------------');
    logmsg(1,'-- Propagation extended to %dth iteration -----------------------', iter);
    % Propagation is being extended, so reset nSteps
    % Continue propagation by 20% more steps or by tcorr/dt, whichever is
    % greater
    Sim.nSteps = max([ceil(tcorrAvg/Sim.dt), ceil(1.2*Sim.nSteps)]);
    Sim.randAngStep = bsxfun(@times, randn(3,Sim.nTraj,Sim.nSteps),...
                                   sqrt(2*Sim.Diff*Sim.dt));
%     Q = propagate(Q, Sim, iter);
    qTraj = propagate(qTraj, Sim, iter);
    
%     % Extend simulation along time axis
%     qTemp = zeros(4,Sim.nTraj,Sim.nSteps);
%     qTemp(1,:,:) = squeeze(real(Q(1,1,:,:)));
%     qTemp(2,:,:) = squeeze(-imag(Q(1,2,:,:)));
%     qTemp(3,:,:) = squeeze(real(Q(2,1,:,:)));
%     qTemp(4,:,:) = squeeze(imag(Q(2,2,:,:)));
%     qTraj = cat(3,qTraj,qTemp);
                                
  end

  if chkcon
    gr = grstat(qTraj);
    converged = all(gr(:)<1.00001);
  else
    converged = 1;
  end

  iter = iter + 1;
  
  if iter>10
    logmsg(1,'Warning: convergence is very slow. Consider increasing\nlength or number of trajectories.')
  end

end

totSteps = size(qTraj,3);

t = linspace(0, totSteps*Sim.dt, totSteps).';

% Convert to rotation matrices
RTraj = quat2rotmat(qTraj);

logmsg(1,'-- Propagation finished --------------------------------------');
logmsg(1,'--------------------------------------------------------------');

%% Final processing
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
  
  case 2  % Output rotation matrices only
    varargout = {t, RTraj};
    
  case 3  % Output rotation matrices and quaternions
    varargout = {t, RTraj, qTraj};
    
end

clear global EasySpinLogLevel

return

%% Helper functions
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
  case 3 % Diff already rhombic
  otherwise
    err = 'Sys.Diff must have 1, 2 or 3 elements (isotropic, axial, rhombic).';
    return
end

end


function q = propagate(q, Sim, iter)
% Propagate quaternions

randAngStep = Sim.randAngStep;
nSteps = Sim.nSteps;
nTraj = Sim.nTraj;
dt = Sim.dt;
Diff = Sim.Diff;
Coefs = Sim.Coefs;
LMK = Sim.LMK;

if iter>0
  % If propagation is being extended, initialize q from the last set
  if ~isempty(Coefs)
    torque = anistorque(LMK,Coefs,q(:,:,end));
%     c2p2 = Coefs(1);
%     c20 = Coefs(1);
%     A = q(1,:,iStep-1) - 1i*q(4,:,iStep-1);
%     B = -q(3,:,iStep-1) - 1i*q(2,:,iStep-1);
%     Ast = q(1,:,iStep-1) + 1i*q(4,:,iStep-1);
%     Bst = -q(3,:,iStep-1) + 1i*q(2,:,iStep-1);
%     X = A.*Bst + Ast.*B;
%     Y = -1i*(A.*Bst - Ast.*B);
%     Z = A.*Ast - B.*Bst;
%     torque = [ -3*c20*Y.*Z;
%                 3*c20*X.*Z;
%                zeros(size(Z))];
%     torque = [-sqrt(6)*c2p2*Y.*Z;
%               -sqrt(6)*c2p2*X.*Z;
%               2*sqrt(6)*c2p2*X.*Y];
    AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,:,1);
  else
    % If there is no orienting potential, then there is no torque to
    % calculate
    AngStep = randAngStep(:,:,1);
  end

  % Calculate size and normalized axis of angular step
  theta = sqrt(sum(AngStep.^2, 1));
  ux = AngStep(1,:,:)./theta;
  uy = AngStep(2,:,:)./theta;
  uz = AngStep(3,:,:)./theta;

  st = sin(theta/2);
  ct = cos(theta/2);

%   U = [    ct, -ux.*st, -uy.*st, -uz.*st; ...
%        ux.*st,      ct,  uz.*st, -uy.*st; ...
%        uy.*st, -uz.*st,      ct,  ux.*st; ...
%        uz.*st,  uy.*st, -ux.*st,      ct];

  % Calculate q for the first time step

  qinit(1,:) =      q(1,:,end).*ct - q(2,:,end).*ux.*st ...
                - q(3,:,end).*uy.*st - q(4,:,end).*uz.*st;
  qinit(2,:) =      q(2,:,end).*ct + q(1,:,end).*ux.*st ...
                - q(4,:,end).*uy.*st + q(3,:,end).*uz.*st;
  qinit(3,:) =      q(3,:,end).*ct + q(4,:,end).*ux.*st ...
                + q(1,:,end).*uy.*st - q(2,:,end).*uz.*st;
  qinit(4,:) =      q(4,:,end).*ct - q(3,:,end).*ux.*st ...
                + q(2,:,end).*uy.*st + q(1,:,end).*uz.*st;

  q = zeros(4,nTraj,nSteps);
  q(:,:,1) = qinit;
end
  
for iStep=2:nSteps
  if ~isempty(Coefs)
    torque = anistorque(LMK,Coefs,q(:,:,iStep-1));
%     c2p2 = Coefs(1);
%     c20 = Coefs(1);
%     A = q(1,:,iStep-1) - 1i*q(4,:,iStep-1);
%     B = -q(3,:,iStep-1) - 1i*q(2,:,iStep-1);
%     Ast = q(1,:,iStep-1) + 1i*q(4,:,iStep-1);
%     Bst = -q(3,:,iStep-1) + 1i*q(2,:,iStep-1);
%     X = A.*Bst + Ast.*B;
%     Y = -1i*(A.*Bst - Ast.*B);
%     Z = A.*Ast - B.*Bst;
%     torque = [ -3*c20*Y.*Z;
%                 3*c20*X.*Z;
%                zeros(size(Z))];
%     torque = [-sqrt(6)*c2p2*Y.*Z;
%               -sqrt(6)*c2p2*X.*Z;
%               2*sqrt(6)*c2p2*X.*Y];
    AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,:,iStep-1);
  else
    % If there is no orienting potential, then there is no torque to
    % calculate
    AngStep = randAngStep(:,:,iStep-1);
  end

  % Calculate size and normalized axis of angular step
  theta = sqrt(sum(AngStep.^2, 1));
  
  ux = AngStep(1,:)./theta;
  uy = AngStep(2,:)./theta;
  uz = AngStep(3,:)./theta;

  st = sin(theta/2);
  ct = cos(theta/2);
  
%   U = [    ct, -ux.*st, -uy.*st, -uz.*st; ...
%        ux.*st,      ct,  uz.*st, -uy.*st; ...
%        uy.*st, -uz.*st,      ct,  ux.*st; ...
%        uz.*st,  uy.*st, -ux.*st,      ct];

  % Perform propagation
  q(1,:,iStep) =      q(1,:,iStep-1).*ct - q(2,:,iStep-1).*ux.*st ...
                - q(3,:,iStep-1).*uy.*st - q(4,:,iStep-1).*uz.*st;
  q(2,:,iStep) =      q(2,:,iStep-1).*ct + q(1,:,iStep-1).*ux.*st ...
                - q(4,:,iStep-1).*uy.*st + q(3,:,iStep-1).*uz.*st;
  q(3,:,iStep) =      q(3,:,iStep-1).*ct + q(4,:,iStep-1).*ux.*st ...
                + q(1,:,iStep-1).*uy.*st - q(2,:,iStep-1).*uz.*st;
  q(4,:,iStep) =      q(4,:,iStep-1).*ct - q(3,:,iStep-1).*ux.*st ...
                + q(2,:,iStep-1).*uy.*st + q(1,:,iStep-1).*uz.*st;
end

end

%% Deprecated

% function Q = propagate(Q, Sim, iter)
% % Propagate quaternion matrix Q
% 
% randAngStep = Sim.randAngStep;
% nSteps = Sim.nSteps;
% nTraj = Sim.nTraj;
% dt = Sim.dt;
% Diff = Sim.Diff;
% Coefs = Sim.Coefs;
% LMK = Sim.LMK;
% 
% if iter>0
%   % If propagation is being extended, initialize Q matrix from the last set
%   if isnumeric(Coefs)
%     torque = anistorque(Q(:,:,:,end), LMK, Coefs);
% %     c2p2 = 2*Coefs(1);
% %     A = Q(1,1,:,iStep-1);
% %     B = Q(1,2,:,iStep-1);
% %     Ast = Q(2,2,:,iStep-1);
% %     Bst = -Q(2,1,:,iStep-1);
% %     X = A.*Bst + Ast.*B;
% %     Y = -1i*(A.*Bst - Ast.*B);
% %     Z = A.*Ast - B.*Bst;
% %     torque = [-sqrt(6)*c2p2*Y.*Z;
% %               -sqrt(6)*c2p2*X.*Z;
% %               2*sqrt(6)*c2p2*X.*Y];
% %     torque = permute(torque,[1,3,2]);
%     AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,:,1);
%   else
%     % If there is no orienting potential, then there is no torque to
%     % calculate
%     AngStep = randAngStep(:,:,1);
%   end
% 
%   % Calculate size and normalized axis of angular step
%   theta = sqrt(sum(AngStep.^2, 1));
%   ux = AngStep(1,:,:)./theta;
%   uy = AngStep(2,:,:)./theta;
%   uz = AngStep(3,:,:)./theta;
% 
%   st = sin(theta/2);
%   ct = cos(theta/2);
% 
%   % Need '1' in the third dimension so that size(U) matches size(Q(:,:,iStep-1,:))
%   U(1,1,:,1) = ct - 1i*uz.*st;
%   U(1,2,:,1) = -(uy + 1i*ux).*st;
%   U(2,1,:,1) = (uy - 1i*ux).*st;
%   U(2,2,:,1) = ct + 1i*uz.*st;
% 
%   % Calculate Q matrices for the first time step
%   Qinit = matmult(Q(:,:,:,end),U(:,:,:,1));
%   Q = zeros(2,2,nTraj,nSteps);
%   Q(:,:,:,1) = Qinit;
% end
%   
% for iStep=2:nSteps
%   if isnumeric(Coefs)
%     torque = anistorque(Q(:,:,:,iStep-1), LMK, Coefs);
% %     c2p2 = 2*Coefs(1);
% %     A = Q(1,1,:,iStep-1);
% %     B = Q(1,2,:,iStep-1);
% %     Ast = Q(2,2,:,iStep-1);
% %     Bst = -Q(2,1,:,iStep-1);
% %     X = A.*Bst + Ast.*B;
% %     Y = -1i*(A.*Bst - Ast.*B);
% %     Z = A.*Ast - B.*Bst;
% %     torque = [-sqrt(6)*c2p2*Y.*Z;
% %               -sqrt(6)*c2p2*X.*Z;
% %               2*sqrt(6)*c2p2*X.*Y];
% %     torque = permute(torque,[1,3,2]);
%     AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,:,iStep-1);
%   else
%     % If there is no orienting potential, then there is no torque to
%     % calculate
%     AngStep = randAngStep(:,:,iStep-1);
%   end
% 
%   % Calculate size and normalized axis of angular step
%   theta = sqrt(sum(AngStep.^2, 1));
%   
%   ux = AngStep(1,:)./theta;
%   uy = AngStep(2,:)./theta;
%   uz = AngStep(3,:)./theta;
% 
%   st = sin(theta/2);
%   ct = cos(theta/2);
% 
%   % Need '1' in the third dimension so that size(U) matches size(Q(:,:,iStep-1,:))
% %   U = zeros(2,2,nTraj);
%   U(1,1,:) = ct - 1i*uz.*st;                  
%   U(1,2,:) = -(uy + 1i*ux).*st;
%   U(2,1,:) = (uy - 1i*ux).*st;
%   U(2,2,:) = ct + 1i*uz.*st;
% 
%   % Perform propagation
%   Q(:,:,:,iStep) = matmult(Q(:,:,:,iStep-1), U);
% end
% 
% end
% 
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