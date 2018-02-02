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

if isfield(Sys,'PseudoPotFun') || isfield(Sys,'ProbDensFun')
  if isfield(Sys,'Coefs')||isfield(Sys,'LMK')
    error('Please choose either PseudoPotFun or Coefs and LMK for an orienting potential.')
  end
  
  if isfield(Sys,'ProbDensFun')
    ProbDensFun = Sys.ProbDensFun;
    idx = ProbDensFun < 1e-14;
    ProbDensFun(idx) = 1e-14;
    PotFun = -log(ProbDensFun); 
  end
  
  if isfield(Sys,'PseudoPotFun'), PotFun = Sys.PseudoPotFun; end
  
%   PotFun = smooth3(PotFun, 'gaussian');
%   PotFun = smoothn(PotFun, 0.5);
  
  alphaGrid = linspace(-pi, pi, size(PotFun,1));
  betaGrid = linspace(0, pi, size(PotFun,2)+2);
  betaGrid = betaGrid(2:end-1);
  gammaGrid = linspace(-pi, pi, size(PotFun,3));

  if any(isnan(PotFun(:)))
    error('At least one NaN detected in log(PseudoPotFun).')
  end
  
  if any(isinf(PotFun(:)))
    error('At least one inf detected in log(PseudoPotFun).')
  end
  
  pidx = [2, 1, 3];
  
  [dx, dy, dz] = gradient_euler(PotFun, alphaGrid, betaGrid, gammaGrid);
  
%   px = smooth3(px, 'gaussian');
%   py = smooth3(py, 'gaussian');
%   pz = smooth3(pz, 'gaussian');
  
  method = 'linear';
  Gradx = griddedInterpolant({alphaGrid, betaGrid, gammaGrid}, dx, method);
  Grady = griddedInterpolant({alphaGrid, betaGrid, gammaGrid}, dy, method);
  Gradz = griddedInterpolant({alphaGrid, betaGrid, gammaGrid}, dz, method);
  Sim.interpGrad = {Gradx, Grady, Gradz};
  
%   clear logPotFun
%   clear px
%   clear py
%   clear pz
else
  Sim.interpGrad = [];
end


% Discrete Monte carlo settings (Par)
% -------------------------------------------------------------------------

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

elseif isfield(Par,'nSteps') && isfield(Par,'tmax')
  % number of steps and max time are given
  tmax = Par.tmax;
  nSteps = Par.nSteps;
  dt = tmax/Sim.nSteps;

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
  if isfield(Sys,'ProbDensFun') || isfield(Sys,'PseudoPotFun')
    [alphaSamples, betaSamples, gammaSamples] = rejectionsample3d(exp(-PotFun), alphaGrid, betaGrid, gammaGrid, Sim.nTraj);
    Omega = [alphaSamples; 
             betaSamples; 
             gammaSamples];
%   elseif isfield(Sys,'Coefs')
  else
    gridPts = linspace(-1, 1, Sim.nTraj);
    gridPhi = zeros(1, Sim.nTraj);
  %   gridPhi = sqrt(pi*Sim.nTraj)*asin(gridPts);
    gridTheta = acos(gridPts);
  %   gridPsi = zeros(1, Sim.nTraj);
    gridPsi = sqrt(pi*Sim.nTraj)*asin(gridPts); % TODO why does this angle, and not Phi, not affect the spectra?
    Omega = [gridPhi; gridTheta; gridPsi];
  end
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
% q0 = euler2quat(Omega,'active');
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
iter = 1;

logmsg(2,'-- Calculating stochastic trajectories -----------------------');

while ~converged
  if iter==1
    %  Pre-calculate angular steps due to random torques
    %  (Eq. 61 from reference, without factor of 1/2)
    Sim = genRandSteps(Sim,Integrator);

    %  Perform stochastic simulations
    %  (Eqs. 48 and 61 in reference)
    qTraj = propagate_classical(qTraj, Sim, iter);

  else
    logmsg(3,'-- Convergence not obtained -------------------------------------');
    logmsg(3,'-- Propagation extended to %dth iteration -----------------------', iter);
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
  
  if iter>15 && converged==0
    logmsg(1,['Warning: restarting trajectory set due to lack of convergence.\n',...
              'Consider increasing length or number of trajectories.\n'])
    iter = 0;
    % re-initialize trajectories
    Sim.nSteps = nSteps;
    q0 = euler2quat(Omega);
    qTraj = zeros(4,Sim.nTraj,Sim.nSteps);
    qTraj(:,:,1) = q0;
  end

end

% clear wignerdquat

totSteps = size(qTraj,3);

t = linspace(0, totSteps*Sim.dt, totSteps).';

% Convert to rotation matrices
RTraj = quat2rotmat(qTraj);

logmsg(2,'-- Propagation finished --------------------------------------');
logmsg(2,'--------------------------------------------------------------');


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

function [dx, dy, dz] = gradient_euler(Data, aGrid, bGrid, gGrid)
% performs the second-order numerical gradient on a 3D dataset as a
% function of Euler angles alpha, beta, gamma

if ~isvector(aGrid) || ~isvector(bGrid) || ~isvector(gGrid)
  error('aGrid, bGrid, and gGrid must be vectors.')
end

if ndims(Data)~=3
  error('Expected a 3-dimensional array for input data.')
end

sizeData = size(Data);

if sizeData(1)~=length(aGrid) || sizeData(2)~=length(bGrid) || sizeData(3)~=length(gGrid)
  error('Dimensions of input data array must match the sizes of the grids.')
end

% create 3D ndgrid from 1D grid vectors
[AGrid, BGrid, GGrid] = ndgrid(aGrid, bGrid, gGrid);

dAlpha = aGrid(2) - aGrid(1);
dBeta = bGrid(2) - bGrid(1);
dGamma = gGrid(2) - gGrid(1);

% calculate second-order numerical derivatives

da = zeros(sizeData);
da(2:end-1,:,:) = Data(3:end,:,:) - Data(1:end-2,:,:);
da(1,:,:) = 4*Data(2,:,:) - 3*Data(1,:,:) - Data(3,:,:);
da(end,:,:) = 3*Data(end,:,:) + Data(end-2,:,:) - 4*Data(end-1,:,:);
da = da/2/dAlpha;

db = zeros(sizeData);
db(:,2:end-1,:) = Data(:,3:end,:) - Data(:,1:end-2,:);
db(:,1,:) = 4*Data(:,2,:) - 3*Data(:,1,:) - Data(:,3,:);
db(:,end,:) = 3*Data(:,end,:) + Data(:,end-2,:) - 4*Data(:,end-1,:);
db = db/2/dBeta;

dg = zeros(sizeData);
dg(:,:,2:end-1) = Data(:,:,3:end) - Data(:,:,1:end-2);
dg(:,:,1) = 4*Data(:,:,2) - 3*Data(:,:,1) - Data(:,:,3);
dg(:,:,end) = 3*Data(:,:,end) + Data(:,:,end-2) - 4*Data(:,:,end-1);
dg = dg/2/dGamma;

% convert to Cartesian

% dx = - csc(BGrid).*cos(GGrid).*da ...
%      + sin(GGrid).*db ...
%      + cot(BGrid).*cos(GGrid).*dg;
% 
% dy = - csc(BGrid).*sin(GGrid).*da ...
%      - cos(GGrid).*db ...
%      + cot(BGrid).*sin(GGrid).*dg;
% 
% dz = dg;

dx = -csc(BGrid).*cos(GGrid).*da ...
    + sin(GGrid).*db ...
    + cot(BGrid).*cos(GGrid).*dg;

dy = csc(BGrid).*sin(GGrid).*da ...
   + cos(GGrid).*db ...
   - cot(BGrid).*sin(GGrid).*dg;

dz = dg;

end

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
interpGrad = Sim.interpGrad;

if ~isempty(Coefs)
  isEigenPot = 1;
else
  isEigenPot = 0;
end

if ~isempty(interpGrad)
  isNumericPot = 1;
else
  isNumericPot = 0;
end

if iter>1
    % If propagation is being extended, initialize q from the last set
  startStep = 1;
else
  % First step has already been initialized by starting orientation, so
  % skip to second step
  startStep = 2;
end

for iStep = startStep:nSteps
  
  if iStep==1
    qLast = q(:,:,end);
    q = zeros(4,nTraj,nSteps);
  else
    qLast = q(:,:,iStep-1);
  end
  
  currRandAngStep = randAngStep(:,:,iStep);
  
  if isEigenPot
    % use Wigner functions of quaternions to calculate torque
    torque = anistorque(LMK, Coefs, qLast);
    AngStep = bsxfun(@times,torque,Diff*dt) + currRandAngStep;
  elseif isNumericPot
    % use orienting pseudopotential functions of Euler angles to calculate
    % torque
    [alpha, beta, gamma] = quat2euler(qLast,'active');
    pxint = interp3fast(interpGrad{1}, alpha, beta, gamma);
    pyint = interp3fast(interpGrad{2}, alpha, beta, gamma);
    pzint = interp3fast(interpGrad{3}, alpha, beta, gamma);
    torque = [-pxint.'; -pyint.'; -pzint.'];
    AngStep = bsxfun(@times,torque,Diff*dt) + currRandAngStep;
  else
    % If there is no orienting potential, then there is no torque to
    % calculate
    AngStep = currRandAngStep;
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
  
  q1 = qLast(1,:);
  q2 = qLast(2,:);
  q3 = qLast(3,:);
  q4 = qLast(4,:);

  q(1,:,iStep) = q1.*ct - q2.*ux.*st - q3.*uy.*st - q4.*uz.*st;
  q(2,:,iStep) = q2.*ct + q1.*ux.*st - q4.*uy.*st + q3.*uz.*st;
  q(3,:,iStep) = q3.*ct + q4.*ux.*st + q1.*uy.*st - q2.*uz.*st;
  q(4,:,iStep) = q4.*ct - q3.*ux.*st + q2.*uy.*st + q1.*uz.*st;
  
%   diff = 1.0-sqrt(sum(q(:,:,iStep).*q(:,:,iStep), 1));
%   
%   if any(abs(diff(:)) > 1e-14)
%     error('Quaternions are not normalized on step %d.\n', iStep)
%   end
  
end

end

function Vq = interp3fast(F, Xq, Yq, Zq)
% Adapted from MATLAB's interp3 function, and makes the following
% assumptions regarding input:
%    Monotonic grid vectors X, Y, Z were fed to griddedInterpolant to
%    obtain F
%    V is ndgrid-ordered, not meshgrid-ordered (fed to griddedInterpolant)
%interp3fast 3-D interpolation (table lookup).

Vq = F([Xq.',Yq.',Zq.']);

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