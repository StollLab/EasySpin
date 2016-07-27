% stochtraj  Generate stochastic rotational trajectories
%
%  [t,R] = stochtraj(Par)
%  [t,R,qtraj] = stochtraj(...)
%
%  Par: structure with simulation parameters
%     .tcorr          correlation time (in seconds, 1 or 3 elements)
%     .lambda         ordering potential coefficient
%     .dt             time step (in seconds)
%     .nSteps         number of time steps per simulation
%     .nTraj          number of trajectories
%     .theta          polar angle of starting orientation (in radians)
%     .phi            azimuthal angle of starting orientation (in radians)
%
%  Output:
%     t              time points of the trajectory (in seconds)
%     R              array of rotation matrices
%     qtraj          array of normalized quaternions

% Implementation based on 
%   Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%     http://dx.doi.org/10.1063/1.2908075

function varargout = stochtraj(Par)

if (nargin == 0), help(mfilename); return; end

if (nargin > 1), error('Too many input arguments.'); end

switch nargout
  case 0 % plotting
  case 2 % t,R
  case 3 % t,R,qTraj
  otherwise
    error('Incorrect number of output arguments.');
end

%========================================================================
% Dynamics and ordering potential
%========================================================================
% if no ordering potential coefficient is given, set it to 0
if ~isfield(Par,'lambda'), Par.lambda = 0; end

if numel(Par.lambda) > 1
  error('Only one orienting potential coefficient, c20, is currently implemented.')
end
lambda = Par.lambda;
if lambda < 0
  error('Orienting potential coefficient(s) cannot be negative.');
end
orderingPotential = (lambda > 0);

% parse the dynamics parameter input using private function
if isfield(Par,'tcorr'), Dynamics.tcorr = Par.tcorr; end
if isfield(Par,'Diff'), Dynamics.Diff = Par.Diff; end
if isfield(Par,'logtcorr'), Dynamics.logtcorr = Par.logtcorr; end
if isfield(Par,'logDiff'), Dynamics.logDiff = Par.logDiff; end

[Dynamics, err] = processdynamics(Dynamics);
error(err);
Diff = Dynamics.Diff';

%========================================================================
% Discrete Monte carlo settings
%========================================================================
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
  t = linspace(0, nSteps*dt, nSteps);

elseif isfield(Par, 'nSteps') && isfield(Par, 'tmax')
  % number of steps and max time are given
  tmax = Par.tmax;
  nSteps = Par.nSteps;
  dt = tmax/nSteps;
  t = linspace(0, nSteps*dt, nSteps);

else
  error(['You must specify a time array, or a number of steps and ' ...
        'either a time step or tmax.'])
end

%========================================================================
% Grid and trajectory settings
%========================================================================
% If number of trajectories is not given, set it to 1
if ~isfield(Par, 'nTraj'), Par.nTraj = 1; end
nTraj = Par.nTraj;

% Get user-supplied starting angles
theta = [];
phi = [];
chi = [];
if isfield(Par, 'theta'), theta = Par.theta; end
if isfield(Par, 'phi'), phi = Par.phi; end
%if isfield(Par, 'chi'), phi = Par.chi; end

% Supplement starting angles if necessary
if isempty(phi)
  phi = rand(1,nTraj)*2*pi;
end
if isempty(theta)
  z = 2*rand(1,nTraj)-1
  theta = acos(z);
end
if isempty(chi)
  % Orienting potentials with m'=0 are independent of chi, so we can set chi0=0
  chi = zeros(1,nTraj);
end

% If only one starting angle and multiple trajectories, repeat the angle
if numel(theta) == 1 && nTraj > 1
  theta = repmat(theta,1,nTraj);
elseif numel(theta) ~= nTraj
  error('The number of starting angles must be equal to the number of trajectories.');
end
if numel(phi) == 1 && nTraj > 1
  phi = repmat(phi,1,nTraj);
elseif numel(phi) ~= nTraj
  error('The number of starting angles must be equal to the number of trajectories.');
end

assert(numel(theta) == numel(phi), 'Theta and phi must be the same size.')
assert(numel(theta) == numel(chi), 'Theta and chi must be the same size.')

q0 = euler2quat(chi,theta,phi);

% Convert initial quaternion to a unitary 2x2 matrix for easier propagation
Q = zeros(2,2,nSteps,nTraj);
Q(1,1,1,:) =  q0(1,:) - 1i*q0(4,:);
Q(1,2,1,:) = -q0(3,:) - 1i*q0(2,:);
Q(2,1,1,:) =  q0(3,:) - 1i*q0(2,:);
Q(2,2,1,:) =  q0(1,:) + 1i*q0(4,:);

% Pre-calculate angular steps due to random torques
%   (Eq. 61 from reference, without factor of 1/2)
randAngStep = bsxfun(@times, randn(3,nSteps-1,nTraj), sqrt(2*Diff*dt));

if ~orderingPotential
  % No ordering potential -> only random torque present
  %  (Eqs. 48 and 61 in reference)
  
  % Calculate size and normalized axis of angular step
  theta = sqrt(sum(randAngStep.^2,1));
  ux = randAngStep(1,:,:)./theta;
  uy = randAngStep(2,:,:)./theta;
  uz = randAngStep(3,:,:)./theta;
  
  % Pre-allocate array of propagators
  U = zeros(2,2,nSteps-1,nTraj);
  
  st = sin(theta/2);
  ct = cos(theta/2);
  U(1,1,:,:) = ct - 1i*uz.*st;
  U(1,2,:,:) = -(uy + 1i*ux).*st;
  U(2,1,:,:) = (uy - 1i*ux).*st;
  U(2,2,:,:) = ct + 1i*uz.*st;
    
  % Perform propagation
  for iStep = 2:nSteps
    Q(:,:,iStep,:) = matmult(Q(:,:,iStep-1,:), U(:,:,iStep-1,:));
  end

else
  % Ordering potential present -> include systematic torque
  %  (Eqs. 48 and 61 in reference)
  for iStep = 2:nSteps
    torque = anistorque(Q(:, :, iStep-1, :), lambda);
    AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,iStep-1,:);

    % Calculate size and normalized axis of angular step
    theta = sqrt(sum(AngStep.^2, 1));
    ux = AngStep(1,:,:)./theta;
    uy = AngStep(2,:,:)./theta;
    uz = AngStep(3,:,:)./theta;
    
    % Pre-allocate and calculate propagator array
    U = zeros(2,2,1,nTraj);
    st = sin(theta/2);
    ct = cos(theta/2);
    U(1,1,:) = ct - 1i*uz.*st;
    U(1,2,:) = -(uy + 1i*ux).*st;
    U(2,1,:) = (uy - 1i*ux).*st;
    U(2,2,:) = ct + 1i*uz.*st;

    % Perform propagation
    Q(:,:,iStep,:) = matmult(Q(:,:,iStep-1,:), U);
  end
  
end

% Convert 2x2 matrices to 4x1 quaternions
qTraj = zeros(4,nSteps,nTraj);
qTraj(1,:,:) = squeeze(real(Q(1,1,:,:)));
qTraj(2,:,:) = squeeze(-imag(Q(1,2,:,:)));
qTraj(3,:,:) = squeeze(real(Q(2,1,:,:)));
qTraj(4,:,:) = squeeze(imag(Q(2,2,:,:)));

% Convert to rotation matrices
R = quat2rotmat(qTraj);

%==============================================================
%  Final processing
%==============================================================

switch nargout
  case 0 % Plot results
    maxTraj = 3;
    if nTraj>maxTraj
      error('Cannot plot more than %d trajectory.',maxTraj);
    end
    clf
    hold on
    for iTraj = 1:min(maxTraj,nTraj)
      x = squeeze(R(1,3,:,iTraj));
      y = squeeze(R(2,3,:,iTraj));
      z = squeeze(R(3,3,:,iTraj));
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
    varargout = {t, R};
    
  case 3  % Output rotation matrices and quaternions
    varargout = {t, R, qTraj};
    
end


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

% if isfield(Dyn,'lw')
%   if numel(Dyn.lw)>1
%     if FieldSweep
%       LorentzFWHM = Dyn.lw(2)*28 * 1e6; % mT -> MHz -> Hz
%     else
%       LorentzFWHM = Dyn.lw(2)*1e6; % MHz -> Hz
%     end
%   else
%     LorentzFWHM = 0;
%   end
%   if (LorentzFWHM~=0)
%     % Lorentzian T2 from FWHM in freq domain 1/T2 = pi*FWHM
%     Dyn.T2 = 1/LorentzFWHM/pi;
%   else
%     Dyn.T2 = inf;
%   end
% end
% 
% % Heisenberg exchange
% %------------------------------------------------------------------
% if ~isfield(Dyn,'Exchange'), Dyn.Exchange = 0; end
% Dyn.Exchange = Dyn.Exchange*2*pi*1e6; % MHz -> angular frequency

return