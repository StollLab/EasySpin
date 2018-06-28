% stochtraj  Generate stochastic rotational trajectories
%
%   [t,qTraj] = stochtraj(Sys)
%   [t,stateTraj] = stochtraj(Sys,Par)
%   [t,qTraj,stateTraj] = stochtraj(...)
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
%     Omega          numeric, size = (3,1) or (3,nTraj)
%                    Euler angles for starting orientation(s) of rotational
%                    diffusion
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
%     Model          string
%                    'Continuous': simulates rotational diffusion with 
%                    continuous degrees of freedom using quaternions
%                    'Discrete': kinetic Monte carlo simulations with 
%                    discrete degrees of freedom (states)
% %                    'Hierarchical': simulates using a hybrid model of both
% %                    discrete (states) and continuous degrees of freedom
% %                    (quaternions)
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
%                    trajectories of states, output from kinetic Monte
%                    carlo simulation

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

if ~isfield(Opt,'Model')
  Model = 'Continuous';
  Opt.statesOnly = 0;
else
  Model = Opt.Model;
  switch Model
    case 'Continuous'
      if all(nargout~=[0,2])
        error('Continous model requires two output arguments.')
      end
      Opt.statesOnly = 0;
    case 'Discrete'
      if ~isfield(Opt,'statesOnly')
        Opt.statesOnly = 0;
      end
      
      if Opt.statesOnly
        if all(nargout~=[0,2])
          error('Discrete model for states only requires two output arguments.')
        end
      else
        if all(nargout~=[0,2,3])
          error('Discrete model with orientations requires two or three output arguments.')
        end
      end
    case 'Hierarchical'
      if ~isfield(Opt,'statesOnly')
        Opt.statesOnly = 0;
      end
      if all(nargout~=[0,2,3])
        error('Hierarchical model requires two or three output arguments.')
      end
    otherwise
      error('Stochastic simulation model not recognized.')
  end
end


global EasySpinLogLevel;
EasySpinLogLevel = Opt.Verbosity;


% Check dynamics and ordering potential
% -------------------------------------------------------------------------

% for stochastic rotational dynamics
if strcmp(Model,'Continuous') || strcmp(Model,'Hierarchical')

  % FieldSweep is not valid for stochtraj, so give empty third arg
  [Dynamics,Sim] = validate_dynord('stochtraj',Sys,[]);

  Sim.Diff = Dynamics.Diff';

  tcorrAvg = 1/6/mean(Dynamics.Diff);

  if isfield(Sys, 'Potential')
      if ~isfield(Sys.Potential, 'ProbDensFun')
        Sys.Potential.ProbDensFun = [];
      end
      if ~isfield(Sys.Potential, 'PseudoPotFun')
        Sys.Potential.PseudoPotFun = [];
      end

    if ~isempty(Sys.Potential.ProbDensFun) || ~isempty(Sys.Potential.PseudoPotFun)
      if isfield(Sys.Potential,'lambda')||isfield(Sys.Potential,'LMK')
        error('Please choose either PseudoPotFun or lambda and LMK for an orienting potential.')
      end
      
      isUserPotFun = 1;

      if ~isempty(Sys.Potential.ProbDensFun)
        ProbDensFun = Sys.Potential.ProbDensFun;
        idx = ProbDensFun < 1e-14;
        ProbDensFun(idx) = 1e-14;
        PseudoPotFun = -log(ProbDensFun); 
      end

      if ~isempty(Sys.Potential.PseudoPotFun), PseudoPotFun = Sys.Potential.PseudoPotFun; end

    %   PotFun = smooth3(PotFun, 'gaussian');
    %   PotFun = smoothn(PotFun, 0.5);

      alphaGrid = linspace(-pi, pi, size(PseudoPotFun,1));
    %   betaGrid = linspace(0, pi, size(PseudoPotFun,2));
      betaGrid = linspace(0, pi, size(PseudoPotFun,2)+2);
      betaGrid = betaGrid(2:end-1);
      gammaGrid = linspace(-pi, pi, size(PseudoPotFun,3));

      if any(isnan(PseudoPotFun(:)))
        error('At least one NaN detected in log(PseudoPotFun).')
      end

      if any(isinf(PseudoPotFun(:)))
        error('At least one inf detected in log(PseudoPotFun).')
      end

      pidx = [2, 1, 3];

      [dx, dy, dz] = gradient_euler(PseudoPotFun, alphaGrid, betaGrid, gammaGrid);

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
      isUserPotFun = 0;
    end
  else
    Sim.interpGrad = [];
    isUserPotFun = 0;
  end
  

else
  Sim.ProbDensFun = [];
  Sim.PseudoPotFun = [];
  isUserPotFun = 0;
end

% for kinetic Monte carlo
if strcmp(Model,'Discrete') || strcmp(Model,'Hierarchical')
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
    error('A transition rate matrix or a transition probability matrix is required for Discrete and Hierarchical Models.')
  end
  
  if ~Opt.statesOnly
    if isfield(Sys,'States')
      States = Sys.States;
      if size(States,1)~=3 || size(States,2)~=nStates
        error(['The size of States must be (3,nStates), with the size of the ' ...
               'second dimension equal to the number of rows (and columns) of Rates.'])
      end
    else
      error('A set of States is required for Discrete and Hierarchical Models for orientations.')
    end
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
if strcmp(Model,'Discrete') || strcmp(Model,'Hierarchical')
  if isfield(Sys,'TransRates')
    TPM = expm(Sim.dt*TRM);
  end
%   TPM = expeig(Sim.dt*Rates);
  cumulTPM = cumsum(TPM,1);
  if any(abs(1-cumulTPM(end,:))>1e-13)
    error('The columns of cumulTPM = sum(TPM,1) must sum to 1.')
  end
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
if strcmp(Model,'Discrete')
  if isfield(Par,'States0')
    States0 = Par.States0;
    if ~all(size(States0)==[1,1])&&~all(size(States0)==[1,Par.nTraj])
      error('States0 should be of size (1,1) or (1,Par.nTraj).')
    end
    if any(States0<1) || any(States0>nStates)
      error('Each entry in States0 needs to be equal to an integer within the range [1,nStates].')
    end
  else
    States0 = randi(nStates,1,Par.nTraj);
  end
end

% Get user-supplied starting angles
Omega = [];
if isfield(Par,'Omega'), Omega = Par.Omega; end

% Supplement starting angles if necessary
if isempty(Omega)
  if isUserPotFun
    [alphaSamples, betaSamples, gammaSamples] = cardamom_rejectionsample3d(exp(-PseudoPotFun), alphaGrid, betaGrid, gammaGrid, Sim.nTraj);
    Omega = [alphaSamples; 
             betaSamples; 
             gammaSamples];
%   elseif isfield(Sys,'lambda')
  elseif strcmp(Model,'Discrete')
    if ~Opt.statesOnly
      Omega = zeros(3,Par.nTraj);
      for iTraj = 1:Par.nTraj
        Omega(:,iTraj) = States(:,States0(iTraj));
      end
    end
  else
    gridPts = linspace(-1, 1, Sim.nTraj);
%     gridPhi = zeros(1, Sim.nTraj);
    gridPhi = sqrt(pi*Sim.nTraj)*asin(gridPts);
    gridTheta = acos(gridPts);
  %   gridPsi = zeros(1, Sim.nTraj);
    gridPsi = sqrt(pi*Sim.nTraj)*asin(gridPts); % TODO why does this angle, and not Phi, not affect the spectra?
    Omega = [gridPhi; gridTheta; gridPsi];
  end
end

% If only one starting angle and multiple trajectories, repeat the angle
if ~Opt.statesOnly
  switch size(Omega,1)
    case 3
      % Euler angles
      if size(Omega,2)==1 && Sim.nTraj > 1
        Omega = repmat(Omega,1,Sim.nTraj);
      elseif size(Omega,2)~=Sim.nTraj
        error('Number of starting orientations must be equal to 1 or nTraj.')
      end
      % set starting quaternions
      % q0 = euler2quat(Omega,'active');
      q0 = euler2quat(Omega);
    case 4
      % quaternions
      if size(Omega,2)==1 && Sim.nTraj > 1
        Omega = repmat(Omega,1,Sim.nTraj);
      elseif size(Omega,2)~=Sim.nTraj
        error('Number of starting orientations must be equal to 1 or nTraj.')
      end
      q0 = Omega;
    otherwise
      error('The size of Omega must be (3,1) or (3,nTraj) for Euler angles, or (4,1) or (4,nTraj) for quaternions.')
  end
end

if isfield(Opt,'chkcon')
  chkcon = Opt.chkcon;
  if chkcon==1 && Sim.nTraj==1
    error('Checking convergence of a single trajectory using the R statistic is not supported.\n')
  end
else
  chkcon = 0;
end

if (strcmp(Model,'Discrete')&&~Opt.statesOnly) || strcmp(Model,'Continuous')
  % initialize quaternion trajectories and their starting orientations
  qTraj = zeros(4,Sim.nTraj,Sim.nSteps);
  qTraj(:,:,1) = q0;
end

if strcmp(Model,'Discrete') && ~Opt.statesOnly
  States = euler2quat(States);
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

switch Model
  case 'Continuous'
    while ~converged
      if iter==1
        %  Perform stochastic simulations
        %  (Eqs. 48 and 61 in reference)
        qTraj = stochtraj_proprottraj(qTraj, Sim, iter);
      else
        logmsg(3,'-- Convergence not obtained -------------------------------------');
        logmsg(3,'-- Propagation extended to %dth iteration -----------------------', iter);
        % Propagation is being extended, so reset nSteps
        % Continue propagation by 20% more steps or by tcorr/dt, whichever is
        % greater
        Sim.nSteps = max([ceil(tcorrAvg/Sim.dt), ceil(1.2*Sim.nSteps)]);
        qTraj = stochtraj_proprottraj(qTraj, Sim, iter);
      end

      if chkcon
        gr = stochtraj_grstat(qTraj);
%         converged = all(gr(:)<1.00001);
        converged = all(gr(:)<1.000001);
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
    totSteps = size(qTraj,3);
  case 'Discrete'
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
  case 'Hierarchical'
    error('Hierarchical model not implemented yet.')
end

    % clear wignerdquat

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
    switch Model
      case 'Continuous'
        varargout = {t, qTraj};
      case 'Discrete'
        if Opt.statesOnly
          varargout = {t, stateTraj};
        else
          varargout = {t, qTraj};
        end
    end
  
  case 3  % Output both quaternion and state trajectories
    switch Model
      case 'Continuous'
        error('Too many output variables for Continuous model.')  % TODO: put this check at the beginning
      case 'Discrete' 
        varargout = {t, qTraj, stateTraj};
      case 'Hierarchical'
        varargout = {t, qTraj, stateTraj};
    end

end

clear global EasySpinLogLevel


end

% Helper functions
% -------------------------------------------------------------------------

function C = expeig(A)

[V,D] = eig(A);

C = V*diag(exp(diag(D)))*V';

end

function [dx, dy, dz] = gradient_euler(Data, aGrid, bGrid, gGrid)
% performs the second-order numerical gradient on a 3D array of Euler 
% angles

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