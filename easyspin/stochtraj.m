% stochtraj  Generate stochastic rotational trajectories
%
%   [t,RTraj] = stochtraj(Sys)
%   [t,RTraj] = stochtraj(Sys,Par)
%   [t,RTraj,qTraj] = stochtraj(...)
%
%   Sys: stucture with system's dynamical parameters
%
%     tcorr          double or numeric, size = (1,3)
%                    correlation time (in seconds)
%
%     logtcorr       double or numeric, size = (1,3)
%                    log10 of rotational correlation time (in seconds)
%
%     Diff           double or numeric, size = (1,3)
%                    diffusion rate (s^-1)
%
%     logDiff        double or numeric, size = (1,3)
%                    log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     Coefs          numeric, size = (nCoefs,2)
%                    array of orienting potential coefficients, with each row
%                    consisting of the corresponding real and imaginary parts
%
%     LMK            numeric, size = (nCoefs,3)
%                    quantum numbers L, M, and K corresponding to each set of 
%                    coefficients
%
%     PseudoPotFun   function handle
%                    orienting pseudopotential function to be used for
%                    calculating the torque
%
%
%   Par: simulation parameters for Monte Carlo integrator
%     dt             double
%                    time step (in seconds)
%
%     nSteps         double
%                    number of time steps per simulation
%
%     nTraj          double
%                    number of trajectories
%
%     Omega          numeric, size = (3,1) or (3,nTraj)
%                    Euler angles for starting orientation(s)
%
%
%   Opt: simulation options
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
%   Output:
%     t              numeric, size = (nSteps,1) 
%                    time points of the trajectory (in seconds)
%
%     RTraj          numeric, size = (3,3,nTraj,nSteps)
%                    rotation matrices
%
%     qTraj          numeric, size = (4,nTraj,nSteps)
%                    normalized quaternions

% Implementation based on 
%  [1] Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%       http://dx.doi.org/10.1063/1.2908075
%  [2] Leimkuhler, Appl.Math.Res.Express 2013, 34 (2013)
%       https://doi.org/10.1093/amrx/abs010
%    newer method for SDE integration that starts with weak convergence of 
%    order 1, but approaches weak convergence of order 2 exponentially fast

function varargout = stochtraj(Sys,Par,Opt)

% Preprocessing
% -------------------------------------------------------------------------

% Only Sys needs to be given to run stochtraj properly, so if Par is not 
% given, initialize it here
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
  case 2 % t,RTraj
  case 3 % t,RTraj,qTraj
  otherwise
    error('Incorrect number of output arguments.');
end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;


% Dynamics and ordering potential
% -------------------------------------------------------------------------

% FieldSweep is not valid for stochtraj, so give empty third arg
[Dynamics,Sim] = validate_dynord('stochtraj',Sys,[]);

Sim.Diff = Dynamics.Diff';

tcorrAvg = 1/6/mean(Dynamics.Diff);

if isfield(Sys,'PseudoPotFun')
  if isfield(Sys,'Coefs')||isfield(Sys,'LMK')
    error('Please choose either PseudoPotFun or Coefs and LMK. Both cannot be used.')
  end
%   if ~isa(PseudoPotFun,'function_handle')
%     error('PseudoPotFun needs to be a function handle.')
%   end
  da = 2*pi/size(Sys.PseudoPotFun,1);
  db = pi/size(Sys.PseudoPotFun,2);
  dg = 2*pi/size(Sys.PseudoPotFun,3);
  
  % Sys.PseudoPotFun may contain zeros, so taking the logarithm directly
  % will yield infs, but since we want to take the gradient of the log of
  % this data anyway, it is better to perform logarithmic differentiation
  % directly and then set infs to zero (which would be due to dividing by 
  % zero, in which case the torque would be negligible at those points 
  % anyway)
  
  idx = Sys.PseudoPotFun<1e-10;
  Sys.PseudoPotFun(idx) = 1e-10;
  
  logpotfun = smooth3(log(Sys.PseudoPotFun), 'gaussian');
  
  [px, py, pz] = gradient(logpotfun, da, db, dg);
  
  Agrid = linspace(-pi, pi, size(Sys.PseudoPotFun,1));
  Bgrid = linspace(0, pi, size(Sys.PseudoPotFun,2));
  Ggrid = linspace(-pi, pi, size(Sys.PseudoPotFun,3));

  extrap = 'none';
  method = 'linear';
  Fx = griddedInterpolant({Agrid, Bgrid, Ggrid}, px, method,extrap);
  Fy = griddedInterpolant({Agrid, Bgrid, Ggrid}, py, method,extrap);
  Fz = griddedInterpolant({Agrid, Bgrid, Ggrid}, pz, method,extrap);
  Sim.interpF = {Fx, Fy, Fz};
else
  Sim.interpF = [];
end


% Discrete Monte carlo settings (Par)
% -------------------------------------------------------------------------

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

elseif isfield(Par,'nSteps') && isfield(Par,'tmax')
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

if ~isfield(Par,'Integrator')
  % default Monte Carlo integrator is Euler-Maruyama
  Integrator = 'Euler-Maruyama';
else
  Integrator = Par.Integrator;
  if ~strcmp(Integrator,'Euler-Maruyama')&&~strcmp(Integrator,'Leimkuhler-Matthews')
    error('Input for integrator method not recognized.')
  end
end



% Grid and trajectory settings
% -------------------------------------------------------------------------

% If number of trajectories is not given, set it to 1
if ~isfield(Par, 'nTraj'), Par.nTraj = 1; end
Sim.nTraj = Par.nTraj;

% Get user-supplied starting angles
Omega = [];
if isfield(Par,'Omega'), Omega = Par.Omega; end

% Supplement starting angles if necessary
if isempty(Omega)
%   if isfield(Sys.PseudoPotFun)
%     
%   else
    z = 2*rand(1,Sim.nTraj)-1;
    Omega = [pi*(2*rand(1,Sim.nTraj)-1);
                            acos(z);
             pi*(2*rand(1,Sim.nTraj)-1)];
%   end
end

% If only one starting angle and multiple trajectories, repeat the angle
if size(Omega,1)==3 
  if size(Omega,2)==1 && Sim.nTraj > 1
    Omega = repmat(Omega,1,Sim.nTraj);
  elseif size(Omega,2)~=Sim.nTraj
    error('Number of starting orientations must be equal to 1 or nTraj.')
  end
else
  error('The size of Omega must be (3,1) or (3,nTraj).')
end

if isfield(Opt,'chkcon')
  chkcon = Opt.chkcon;
  if chkcon==1 && Sim.nTraj==1
    error('Checking convergence of a single trajectory using the R statistic is not supported.\n')
  end
else
  chkcon = 0;
end

% initialize quaternion trajectories and their starting orientations
q0 = euler2quat(Omega);
qTraj = zeros(4,Sim.nTraj,Sim.nSteps);
qTraj(:,:,1) = q0;

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
iter = 0;

logmsg(1,'-- Calculating stochastic trajectories -----------------------');

while ~converged
  if iter==0
    %  Pre-calculate angular steps due to random torques
    %  (Eq. 61 from reference, without factor of 1/2)
    Sim = genRandSteps(Sim,Integrator);

    %  Perform stochastic simulations
    %  (Eqs. 48 and 61 in reference)
    qTraj = propagate_classical(qTraj, Sim, iter);

  else
    logmsg(1,'-- Convergence not obtained -------------------------------------');
    logmsg(1,'-- Propagation extended to %dth iteration -----------------------', iter);
    % Propagation is being extended, so reset nSteps
    % Continue propagation by 20% more steps or by tcorr/dt, whichever is
    % greater
    Sim.nSteps = max([ceil(tcorrAvg/Sim.dt), ceil(1.2*Sim.nSteps)]);
    Sim = genRandSteps(Sim,Integrator);
    qTraj = propagate_classical(qTraj, Sim, iter);

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

% clear wignerdquat

totSteps = size(qTraj,3);

t = linspace(0, totSteps*Sim.dt, totSteps).';

% Convert to rotation matrices
RTraj = quat2rotmat(qTraj);

logmsg(1,'-- Propagation finished --------------------------------------');
logmsg(1,'--------------------------------------------------------------');


% Final processing
% -------------------------------------------------------------------------

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


end

% Helper functions
% -------------------------------------------------------------------------

function Sim = genRandSteps(Sim,Integrator)
% generate random angular steps using one of two different MC integrators

% generate Gaussian random deviates
randns = randn(3,Sim.nTraj,Sim.nSteps);

if strcmp(Integrator,'Euler-Maruyama')
  % standard method for integrating a first order SDE
  Sim.randAngStep = bsxfun(@times, randns, sqrt(2*Sim.Diff*Sim.dt));
elseif strcmp(Integrator,'Leimkuhler-Matthews')
  % see Ref. [2]
  Sim.randAngStep = bsxfun(@times, (randns(:,:,1:end-1)+randns(:,:,2:end))/2, ...
                                   sqrt(2*Sim.Diff*Sim.dt));
end

end

function q = propagate_classical(q, Sim, iter)
% Propagate quaternions

randAngStep = Sim.randAngStep;
nSteps = Sim.nSteps;
nTraj = Sim.nTraj;
dt = Sim.dt;
Diff = Sim.Diff;
Coefs = Sim.Coefs;
LMK = Sim.LMK;

if ~isempty(Sim.interpF)
  interpF = Sim.interpF;
end

if iter>0
  % If propagation is being extended, initialize q from the last set
  if ~isempty(Coefs)
    % use Wigner functions of quaternions to calculate torque
    torque = anistorque(LMK, Coefs, q(:,:,end));
    AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,:,1);
  elseif ~isempty(interpF)
    % use orienting pseudopotential functions of Euler angles to calculate
    % torque
    [alpha, beta, gamma] = quat2euler(q(:,:,end));
    pxint = interp3fast(interpF{1}, alpha, beta, gamma);
    pyint = interp3fast(interpF{2}, alpha, beta, gamma);
    pzint = interp3fast(interpF{3}, alpha, beta, gamma);
    torque = [pxint; pyint; pzint];
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

  qinit(1,:,1) =      q(1,:,end).*ct - q(2,:,end).*ux.*st ...
                - q(3,:,end).*uy.*st - q(4,:,end).*uz.*st;
  qinit(2,:,1) =      q(2,:,end).*ct + q(1,:,end).*ux.*st ...
                - q(4,:,end).*uy.*st + q(3,:,end).*uz.*st;
  qinit(3,:,1) =      q(3,:,end).*ct + q(4,:,end).*ux.*st ...
                + q(1,:,end).*uy.*st - q(2,:,end).*uz.*st;
  qinit(4,:,1) =      q(4,:,end).*ct - q(3,:,end).*ux.*st ...
                + q(2,:,end).*uy.*st + q(1,:,end).*uz.*st;

  q = zeros(4,nTraj,nSteps);
  q(:,:,1) = qinit;
end
  
for iStep=2:nSteps
  if ~isempty(Coefs)
    % use Wigner functions of quaternions to calculate torque
    torque = anistorque(LMK, Coefs, q(:,:,iStep-1));
    AngStep = bsxfun(@times,torque,Diff*dt) + randAngStep(:,:,iStep-1);
  elseif ~isempty(interpF)
    % use orienting pseudopotential functions of Euler angles to calculate
    % torque
    [alpha, beta, gamma] = quat2euler(q(:,:,iStep-1));
    pxint = interp3fast(interpF{1}, alpha, beta, gamma);
    pyint = interp3fast(interpF{2}, alpha, beta, gamma);
    pzint = interp3fast(interpF{3}, alpha, beta, gamma);
    torque = [pxint.'; pyint.'; pzint.'];
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

%   Perform propagation
  q(1,:,iStep) =      q(1,:,iStep-1).*ct - q(2,:,iStep-1).*ux.*st ...
                - q(3,:,iStep-1).*uy.*st - q(4,:,iStep-1).*uz.*st;
  q(2,:,iStep) =      q(2,:,iStep-1).*ct + q(1,:,iStep-1).*ux.*st ...
                - q(4,:,iStep-1).*uy.*st + q(3,:,iStep-1).*uz.*st;
  q(3,:,iStep) =      q(3,:,iStep-1).*ct + q(4,:,iStep-1).*ux.*st ...
                + q(1,:,iStep-1).*uy.*st - q(2,:,iStep-1).*uz.*st;
  q(4,:,iStep) =      q(4,:,iStep-1).*ct - q(3,:,iStep-1).*ux.*st ...
                + q(2,:,iStep-1).*uy.*st + q(1,:,iStep-1).*uz.*st;
  
%   qnorm = sum(q.*q,1);
%   if any(qnorm(:)-1>1e-13)
%     error('qs are not normalized on step %d! Max norm: %e \n', iStep, max(qnorm(:)-1))
%   end

end

end

function Vq = interp3fast(F, Xq, Yq, Zq)
% Adapted from MATLAB's interp3 function, and makes the following
% assumptions regarding input:
%    Output from griddedInterpolant is fed into F
%    Monotonic grid vectors X, Y, Z (fed to griddedInterpolant)
%    V is ndgrid-ordered, not meshgrid-ordered (fed to griddedInterpolant)
%interp3fast 3-D interpolation (table lookup).

p = [2 1 3];

Xq = permute(Xq,p);
Yq = permute(Yq,p);
Zq = permute(Zq,p);

Vq = F(Xq,Yq,Zq);

% we don't want meshgrid-order!
% if iscompact ||  transposedquery
%     % Compact grid evaluation produces a NDGRID
%     % Convert to MESHGRID
%     Vq = permute(Vq,p);
% end

end

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