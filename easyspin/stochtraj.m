% stochtraj  Generate stochastic trajectories of quaternions.
%
%  R = stochtraj(tcorr, lambda, dt, nSteps, nTraj, theta, phi);
%  [R, qtraj] = stochtraj(...)
%
%  Par: simulation parameters
%     tcorr          correlation time (in seconds, 1 or 3 elements)
%     lambda         ordering potential coefficient
%                    Only c20 is currently implemented!
%     dt             time step (in seconds)
%     nSteps         number of time steps per simulation
%     nSims          number of simulations
%     theta          polar angle of starting orientation (in radians)
%     phi            azimuthal angle of starting orientation (in radians)
%
%  Output:
%     R              array of rotation matrices
%     qtraj          array of normalized quaternions
%     t              time points of the trajectory

% Implementation based on 
% - Sezer, et al., J.Chem.Phys. 128, 165106 (2008), doi: 10.1063/1.2908075

function varargout = stochtraj(Par)

if (nargin == 0), help(mfilename); return; end

if (nargin > 1), error('Too many input arguments.'); end
if (nargout > 3), error('Too many output arguments.'); end

%========================================================================
% Dynamics and ordering potential
% Borrowed from chili
%========================================================================
% if no ordering potential coefficient is given, set it to 0
if ~isfield(Par, 'lambda'), Par.lambda = 0; end

if numel(Par.lambda) > 1
  error('Only one orienting potential coefficient, c20, is currently implemented.')
end
lambda = Par.lambda;

% parse the dynamics parameter input using private function
if isfield(Par, 'tcorr'), Dynamics.tcorr = Par.tcorr; end
if isfield(Par, 'Diff'), Dynamics.Diff = Par.Diff; end
if isfield(Par, 'logtcorr'), Dynamics.logtcorr = Par.logtcorr; end
if isfield(Par, 'logDiff'), Dynamics.logDiff = Par.logDiff; end

[Dynamics, err] = processdynamics(Dynamics);
error(err);
Diff = Dynamics.Diff';

%========================================================================
% Discrete Monte carlo settings
%========================================================================
% time axis is given explicitly
if isfield(Par, 't')
  t = Par.t;
  nSteps = numel(t);
  dt = t(2) - t(1);
  if (abs(dt - max(t)/nSteps) > eps), error('t is not linearly spaced.'); end
% number of steps and time step are given
elseif isfield(Par, 'nSteps') && isfield(Par,'dt')
  dt = Par.dt;
  nSteps = Par.nSteps;
  t = linspace(0, nSteps*dt, nSteps);
% number of steps and max time are given
elseif isfield(Par, 'nSteps') && isfield(Par, 'tmax')
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
% if number of trajectories is not given, set it to 1
if ~isfield(Par, 'nTraj'), Par.nTraj = 1; end
if isfield(Par, 'theta') && isfield(Par, 'phi') % && isfield(Par, 'chi')
  theta = Par.theta;
  phi = Par.phi;
  nTraj = Par.nTraj;
  assert(numel(theta) == numel(phi), 'Theta and phi must be the same size.')
  % if only one starting angle and multiple trajectories, repeat the angle
  if numel(theta) == 1 && nTraj > 1
    theta = theta*ones(1, nTraj);
    phi = phi*ones(1, nTraj);
  elseif numel(theta) ~= nTraj
    error(['If the number of starting angles is greater than one, ' ...
           'then it must be equal to the number of trajectories.'])
  end
% if no starting angles are given, generate them from a uniformly random
% distribution
else
  nTraj = Par.nTraj;
  [phi, theta] = sphrand(nTraj);
end

% Note that orienting potentials with m'=0 are independent of alpha, so
% WLOG we can set alpha0=0
% qinit = reshape(euler2quat(zeros(1, numel(theta)), theta, phi), ...
%                 4, 1, nTraj);

qTraj = zeros(4, nSteps, nTraj);
qTraj(:, 1, :) = reshape(euler2quat(zeros(1, numel(theta)), theta, phi), ...
                4, 1, nTraj);

% Convert initial quaternion to a unitary 2x2 matrix for easier propagation
Q = zeros(2, 2, nSteps, nTraj);
Q(1, 1, 1, :) = qTraj(1, 1, :)-1i*qTraj(4, 1, :);
Q(1, 2, 1, :) = -qTraj(3, 1, :)-1i*qTraj(2, 1, :);
Q(2, 1, 1, :) = qTraj(3, 1, :)-1i*qTraj(2, 1, :);
Q(2, 2, 1, :) = qTraj(1, 1, :)+1i*qTraj(4, 1, :);
% Q(:, :, 1, :) = repmat([qinit(1)-1i*qinit(4), -qinit(3)-1i*qinit(2);...
%                         qinit(3)-1i*qinit(2),  qinit(1)+1i*qinit(4)],...
%                         1, 1, 1, nTraj);

% Pre-allocate angular steps due to random torques
% Eq. 61 from reference, without factor of 1/2
randAngStep = bsxfun(@times, randn(3, nSteps-1, nTraj), sqrt(2*Diff*dt));

% randTorque = randn(3, nSteps-1, nSims) ...
%               .*repmat(sqrt(0.5*DiffTensor*dt), 1, nSteps-1, nSims);

if lambda == 0
  % No need to calculate torque due to an ordering potential
  % Eqs. 48 and 61 in reference
  theta = sqrt(sum(randAngStep.^2, 1));
  ux = randAngStep(1, :, :)./theta;
  uy = randAngStep(2, :, :)./theta;
  uz = randAngStep(3, :, :)./theta;
  
  % Pre-allocate propagator
  U = zeros(2, 2, nSteps-1, nTraj);
  
  st = sin(theta/2);
  ct = cos(theta/2);
  U(1, 1, :, :) = ct - 1i*uz.*st;
  U(1, 2, :, :) = -(uy + 1i*ux).*st;
  U(2, 1, :, :) = (uy - 1i*ux).*st;
  U(2, 2, :, :) = ct + 1i*uz.*st;
    
  % Perform propagation and update quaternion trajectory
  for iStep = 2:nSteps
    Q(:, :, iStep, :) = matmult(Q(:, :, iStep-1, :), ...
                                U(:, :, iStep-1, :));
  end

elseif lambda > 0
  % Eqs. 48 and 61 in reference
  for iStep = 2:nSteps
    torque = anistorque(Q(:, :, iStep-1, :), lambda);
    AngStep = bsxfun(@times, torque, Diff*dt) ...
              + randAngStep(:, iStep-1, :);

    % Calculate size and normalized axis of angular step
    theta = sqrt(sum(AngStep.^2, 1));
    ux = AngStep(1, :, :)./theta;
    uy = AngStep(2, :, :)./theta;
    uz = AngStep(3, :, :)./theta;

    % Calculate propagator
    U = zeros(2, 2, 1, nTraj);

    st = sin(theta/2);
    ct = cos(theta/2);
    U(1, 1, :) = ct - 1i*uz.*st;
    U(1, 2, :) = -(uy + 1i*ux).*st;
    U(2, 1, :) = (uy - 1i*ux).*st;
    U(2, 2, :) = ct + 1i*uz.*st;

    % Perform propagation and update quaternion trajectory
    Q(:, :, iStep, :) = matmult(Q(:, :, iStep-1, :), U);
  end
else
  error('Orienting potential coefficient(s) must be positive.')
end

qTraj(1, :, :) = squeeze(real(Q(1, 1, :, :)));
qTraj(2, :, :) = squeeze(-imag(Q(1, 2, :, :)));
qTraj(3, :, :) = squeeze(real(Q(2, 1, :, :)));
qTraj(4, :, :) = squeeze(imag(Q(2, 2, :, :)));

R = quat2rotmat(qTraj);

%==============================================================
%  Final processing
%==============================================================

switch nargout
    
  case 1  % Only output rotation matrices
    varargout = {t, R};
        
  case 2  % Output rotation matrices and quaternions
    varargout = {t, R, qTraj};

end


return

%========================================================================
%========================================================================
% Borrowed from chili for consistency
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