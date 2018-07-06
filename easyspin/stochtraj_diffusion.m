% stochtraj_diffusion  Generate stochastic rotational trajectories
%
%   [t,RTraj] = stochtraj_diffusion(Sys)
%   [t,RTraj,qTraj] = stochtraj_diffusion(...)
%   ... = stochtraj_diffusion(Sys,Par)
%   ... = stochtraj_diffusion(Sys,Par,Opt)
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
%     Potential      structure 
%                    defines an orienting potential using the following 
%                    fields:
%
%       lambda         numeric, size = (nCoefs,2)
%                      array of orienting potential coefficients, with each 
%                      row consisting of the corresponding real and 
%                      imaginary parts
%
%       LMK            numeric, size = (nCoefs,3)
%                      quantum numbers L, M, and K corresponding to each 
%                      set of coefficients
%
%       ProbDensFun    numeric, 3D array
%                      probability distribution grid to be used for
%                      calculating the pseudopotential and the torque
%
%       PseudoPotFun   numeric, 3D array
%                      orienting pseudopotential grid to be used for
%                      calculating the torque
%
%
%
%   Par: simulation parameters for Monte Carlo integrator
%
%     nTraj          integer
%                    number of trajectories
%
%     Omega          numeric, size = (3,1) or (3,nTraj)
%                    Euler angles for starting orientation(s) of rotational
%                    diffusion
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
%     checkConvergence  if set to true, after the first nSteps of the 
%                       trajectories are calculated, both inter- and intra-
%                       trajectory convergence is checked using the Gelman-
%                       Rubin R statistic such that R<1.1, and if this 
%                       condition is not satisfied, then propagation will be 
%                       extended by either a length of time equal to the 
%                       average of tcorr or by 20% more time steps, whichever 
%                       is larger
%
%     Verbosity         0: no display, 1: show info
%
%
%   Output:
%
%     t              matrix, size = (nSteps,1) 
%                    time points of the trajectory (in seconds)
%
%     RTraj          3D array, size = (4,nTraj,nSteps)
%                    trajectories of rotation matrices
%
%     qTraj          3D array, size = (4,nTraj,nSteps)
%                    trajectories of normalized quaternions

function varargout = stochtraj_diffusion(Sys,Par,Opt)

% Preprocessing
%-------------------------------------------------------------------------------

% Only Sys needs to be given to run stochtraj properly, so if Par is not 
% given, initialize it here
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
  case 2 % t, RTraj
  case 3 % t, RTraj, qTraj
  otherwise
    error('Incorrect number of output arguments.');
end

% Check Opt
%-------------------------------------------------------------------------------

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;


% Check dynamics and ordering potential
%-------------------------------------------------------------------------------

% FieldSweep is not valid for stochtraj, so give empty third arg
[Dynamics,Sim] = validate_dynord('stochtraj_diffusion',Sys);
Sim.Diff = Dynamics.Diff;

tcorrAvg = 1/6/mean(Dynamics.Diff);

% Supplement fields
if ~isfield(Sys,'Potential'), Sys.Potential = struct; end
if ~isfield(Sys.Potential,'ProbDensFun'), Sys.Potential.ProbDensFun = []; end
if ~isfield(Sys.Potential,'PseudoPotFun'), Sys.Potential.PseudoPotFun = []; end
if ~isfield(Sys.Potential,'LMK'), Sys.Potential.LMK = []; end
if ~isfield(Sys.Potential,'lambda'), Sys.Potential.lambda = []; end

isUserPotFun = ~isempty(Sys.Potential.ProbDensFun) || ~isempty(Sys.Potential.PseudoPotFun);
isLMKexpansion = ~isempty(Sys.Potential.lambda) || ~isempty(Sys.Potential.LMK);

if isUserPotFun && isLMKexpansion
  error('Please choose either PseudoPotFun or lambda and LMK for an orienting potential.')
end

if isUserPotFun
  
  % Get potential function
  if ~isempty(Sys.Potential.PseudoPotFun)
    PseudoPotFun = Sys.Potential.PseudoPotFun;
  elseif ~isempty(Sys.Potential.ProbDensFun)
    ProbDensFun = Sys.Potential.ProbDensFun;
    ProbDensFun(ProbDensFun < 1e-14) = 1e-14;
    PseudoPotFun = -log(ProbDensFun);
  end
  
  % Check potential function
  if any(isnan(PseudoPotFun(:)))
    error('NaN detected in PseudoPotFun.');
  end  
  if any(isinf(PseudoPotFun(:)))
    error('At least one inf detected in PseudoPotFun.');
  end

  % Calculate gradient
  alphaGrid = linspace(-pi,pi,size(PseudoPotFun,1));
  betaGrid = linspace(0,pi,size(PseudoPotFun,2)+2);
  betaGrid = betaGrid(2:end-1); % avoid poles
  gammaGrid = linspace(-pi,pi,size(PseudoPotFun,3));
  [dx,dy,dz] = gradient_euler(PseudoPotFun, alphaGrid, betaGrid, gammaGrid);
  
  % Set  up interpolant for gradient
  method = 'linear';
  Gradx = griddedInterpolant({alphaGrid, betaGrid, gammaGrid}, dx, method);
  Grady = griddedInterpolant({alphaGrid, betaGrid, gammaGrid}, dy, method);
  Gradz = griddedInterpolant({alphaGrid, betaGrid, gammaGrid}, dz, method);
  Sim.interpGrad = {Gradx, Grady, Gradz};
  
end


% Discrete Monte Carlo settings (Par)
%-------------------------------------------------------------------------

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

elseif isfield(Par,'nSteps') && isfield(Par,'tMax')
  % number of steps and max time are given
  tMax = Par.tMax;
  nSteps = Par.nSteps;
  dt = tMax/Sim.nSteps;

else
  logmsg(1,'-- No time step given. Par.dt set to Par.tcorr/10: %0.5g s.', tcorrAvg);
  dt = tcorrAvg/10;
  if ~isfield(Par,'nSteps')
    nSteps = ceil(200e-9/dt);
    logmsg(1,'-- Number of time steps not given. Par.nSteps set to 200e-9/Par.dt: %d.', nSteps);
  else
    nSteps = Par.nSteps;
  end
end

Sim.nSteps = nSteps;
Sim.dt = dt;

% Set integrator type
if ~isfield(Par,'Integrator')
  % default Monte Carlo integrator is Euler-Maruyama
  Par.Integrator = 'Euler-Maruyama';
end
if ~strcmp(Par.Integrator,'Euler-Maruyama') && ~strcmp(Par.Integrator,'Leimkuhler-Matthews')
  error('Input for integrator method not recognized.')
end
Sim.Integrator = Par.Integrator;


% Grid and trajectory settings
%-------------------------------------------------------------------------------

% If number of trajectories is not given, set it to 1
if ~isfield(Par,'nTraj'), Par.nTraj = 1; end
Sim.nTraj = Par.nTraj;

% Get user-supplied starting angles
Omega = [];
if isfield(Par,'Omega'), Omega = Par.Omega; end

% Supplement starting angles if necessary
if isempty(Omega)
  if isUserPotFun
    % TODO: also do this if lambda+LMK are given
    [alphaSamples, betaSamples, gammaSamples] = cardamom_rejectionsample3d(exp(-PseudoPotFun), alphaGrid, betaGrid, gammaGrid, Sim.nTraj);
    Omega = [alphaSamples; betaSamples; gammaSamples];
  else
    % Set up spiral grid over
    gridPts = linspace(-1, 1, Sim.nTraj);
    gridTheta = acos(gridPts);
    gridPhi = sqrt(pi*Sim.nTraj)*asin(gridPts);
    gridPsi = sqrt(pi*Sim.nTraj)*asin(gridPts); % TODO: why does this angle, and not Phi, not affect the spectra?
    Omega = [gridPhi; gridTheta; gridPsi];
  end
end

% Assure it's a column vector if a vector is given
if isvector(Omega), Omega = Omega(:); end

% If only one starting angle and multiple trajectories, repeat the angle
if size(Omega,2)==1 && Sim.nTraj>1
  Omega = repmat(Omega,1,Sim.nTraj);
end  
  
if size(Omega,2)~=Sim.nTraj
  error('Number of starting orientations must be equal to Par.nTraj.')
end

switch size(Omega,1)
  case 3 % Euler angles
    q0 = euler2quat(Omega);
  case 4 % quaternions
    q0 = Omega;
  otherwise
    error('The size of Omega must be (3,1) or (3,nTraj) for Euler angles, or (4,1) or (4,nTraj) for quaternions.')
end

% initialize quaternion trajectories and their starting orientations
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

if ~isfield(Opt,'convTolerance')
  Opt.convTolerance = 1e-6;
end
if isfield(Opt,'checkConvergence')
  checkConvergence = Opt.checkConvergence;
  if checkConvergence && Sim.nTraj==1
    error('Checking convergence of a single trajectory using the R statistic is not supported.');
  end
else
  checkConvergence = false;
end


% Simulation
%-------------------------------------------------------------------------------

converged = false;

% If checkConvergence=true, then we need to keep track of how many iterations
% have been completed (i.e. the number of instances of propagation in time by nSteps)
iter = 1;

logmsg(2,'-- Calculating stochastic trajectories -----------------------');

while ~converged
  if iter==1
    %  Perform stochastic simulations (Eqs. 48 and 61 in reference)
    qTraj = stochtraj_proprottraj(qTraj, Sim, iter);
  else
    logmsg(3,'-- Convergence not obtained -------------------------------------');
    logmsg(3,'-- Propagation extended to %dth iteration -----------------------', iter);
    % Propagation is being extended, so reset nSteps
    % Continue propagation by 20% more steps or by tcorr/dt, whichever is greater
    Sim.nSteps = max([ceil(tcorrAvg/Sim.dt), ceil(1.2*Sim.nSteps)]);
    qTraj = stochtraj_proprottraj(qTraj, Sim, iter);
  end

  if checkConvergence
    gr = stochtraj_grstat(qTraj);
    converged = all(gr(:) < 1 + Opt.convTolerance);
  else
    converged = true;
  end

  iter = iter + 1;

  if iter>15 && ~converged
    logmsg(1,['Warning: restarting trajectory set due to lack of convergence.\n',...
              'Consider increasing length or number of trajectories.\n'])
    iter = 1;
    % re-initialize trajectories
    Sim.nSteps = nSteps;
    q0 = euler2quat(Omega);
    qTraj = zeros(4,Sim.nTraj,Sim.nSteps);
    qTraj(:,:,1) = q0;
  end
  
end
totSteps = size(qTraj,3);

t = linspace(0, totSteps*Sim.dt, totSteps).';
RTraj = quat2rotmat(qTraj);

logmsg(2,'-- Propagation finished --------------------------------------');
logmsg(2,'--------------------------------------------------------------');


% Final processing
%-------------------------------------------------------------------------------

switch nargout
  case 0 % Plot results
    maxTraj = 3;
    if Sim.nTraj>maxTraj
      error('Cannot plot more than %d trajectories.',maxTraj);
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

  case 2  % Output rotation matrix trajectories
    varargout = {t, RTraj};

  case 3  % Output rotation matrix and quaternion trajectories
    varargout = {t, RTraj, qTraj};
    
end

clear global EasySpinLogLevel


end

% Helper functions
%-------------------------------------------------------------------------------

function [dx, dy, dz] = gradient_euler(Data, aGrid, bGrid, gGrid)
% performs the second-order numerical gradient on a 3D array of Euler angles

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
[~, BGrid, GGrid] = ndgrid(aGrid, bGrid, gGrid);

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
dx = -csc(BGrid).*cos(GGrid).*da ...
    + sin(GGrid).*db ...
    + cot(BGrid).*cos(GGrid).*dg;

dy = csc(BGrid).*sin(GGrid).*da ...
   + cos(GGrid).*db ...
   - cot(BGrid).*sin(GGrid).*dg;

dz = dg;

end

%-------------------------------------------------------------------------------

%{
function varargout = acorr_convergence(RTraj, tol)
% Calculate angular histogram of trajectories and compare with analytic
% expression

nBins = 50;  % There might be a better value to choose here
  
VecTraj = squeeze(RTraj(:, 3, :, :));

totTraj = size(RTraj,4);

AutoCorrFFT = zeros(nSteps, nTraj);

for k = 1:totTraj
  AutoCorrFFT(:, k) = autocorrfft(VecTraj(:, :, k).^2);
end

AutoCorrFFT = sum(AutoCorrFFT, 2)'/totTraj;

converged = ChiSquare < tol;

varargout = {converged, ChiSquare};

end
%}

% -------------------------------------------------------------------------

%{
function varargout = hist_convergence(RTraj, lambda, tol)
% Calculate angular histogram of trajectories and compare with analytic
% expression

nBins = 50;  % There might be a better value to choose here
  
VecTraj = squeeze(RTraj(:, 3, :, :));

totTraj = size(RTraj,4);

bins = linspace(0, pi, nBins)';
ThetaHist = zeros(nBins, totTraj);

for k = 1:totTraj
  ThetaHist(:, k) = hist(acos(VecTraj(3, :, k)), bins);
end

ThetaHist = sum(ThetaHist, 2);
ThetaHist = ThetaHist/sum(ThetaHist);

BoltzDist = exp(lambda*(1.5*cos(bins).^2 - 0.5));
BoltzInt = sum(BoltzDist.*sin(bins));
BoltzDist = BoltzDist.*sin(bins)./BoltzInt;

ChiSquare = sum(((ThetaHist - BoltzDist).^2)./ThetaHist);

converged = ChiSquare < tol;

varargout = {converged, ChiSquare};

end
%}