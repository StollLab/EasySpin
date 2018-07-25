% cardamom  Trajectory-based simulation of CW-EPR spectra.
%
%   cardamom(Sys,Exp,Par)
%   cardamom(Sys,Exp,Par,Opt)
%   cardamom(Sys,Exp,Par,Opt,MD)
%   spc = cardamom(...)
%   [B,spc] = cardamom(...)
%   [B,spc,ExpectVal,t] = cardamom(...)
%
%   Computes a CW-EPR spectrum of an S=1/2 spin label using stochastic or
%   molecular dynamics trajectories.
%
%   Sys: stucture with system's dynamical parameters
%
%     tcorr          double or numeric vector, size = (1,3)
%                    correlation time (in seconds)
%
%     logtcorr       double or numeric vector, size = (1,3)
%                    log10 of rotational correlation time (in seconds)
%
%     Diff           double or numeric vector, size = (1,3)
%                    diffusion rate (s^-1)
%
%     logDiff        double or numeric vector, size = (1,3)
%                    log10 of diffusion rate (s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     DiffGlobal     double or numeric vector, size = (1,3)
%                    global diffusion rate (s^-1)
%
%     Potential      defines an orienting potential using one of the following:
%
%                    a) set of LMKs and lambdas for an expansion in terms of
%                       Wigner functions [L1 M1 K1 lambda1; L2 M2 K2 lambda2]
%                       lambda can be real or complex
%                    b) 3D array defining the potential over a grid, with the
%                       first dimension representing alpha [0,2*pi], the second
%                       beta [0,pi], and the third gamma [0,2*pi]
%                    c) a function handle for a function that takes three
%                       arguments (alpha, beta, and gamma) and returns the value
%                       of the ordering potential for that orientation. The
%                       function should be vectorized, i.e. work with arrays of
%                       alpha, beta, and gamma.
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
%     Sys.lw         double or numeric vector, size = (1,2)
%                    vector with FWHM residual broadenings
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
%
%
%   Par: structure with simulation parameters
%     dt             double
%                    rotational dynamics propagation time step (in seconds)
%
%     Dt             double
%                    spin dynamics propagation time step (in seconds)
%
%     nSteps         integer
%                    number of time steps per simulation
%
%     nTraj          integer
%                    number of trajectories
%
%     OriStart       numeric, size = (3,1), (1,3), or (3,nTraj)
%                    Euler angles for starting orientation(s)
%
%     nOrients       integer
%                    number of lab-to-molecule orientations to loop over
%
%     Orients        numeric matrix, size = (2,nOrients)
%                    (optional) (phi,theta) angles of lab-to-molecule 
%                    orientations, if not given, these are chosen as points
%                    on a spherical spiral grid
%
%     Model          string: 'stochastic', 'jump', 'MD'
%
%
%   Exp: experimental parameter settings
%     mwFreq         double
%                    microwave frequency, in GHz (for field sweeps)
%
% %     Range          numeric vector, size = (1,2)
% %                    sweep range, [sweepmin sweepmax], in mT (for field sweep)
% %
% %     CenterSweep    numeric vector, size = (1,2)
% %                    sweep range, [center sweep], in mT (for field sweep)
% % 
% %     nPoints        integer
% %                    number of points
% % 
% %     Harmonic       integer
% %                    detection harmonic: 0, 1 (default), 2
% % 
% %     ModAmp         double
% %                    peak-to-peak modulation amplitude, in mT (field sweeps only)
% % 
% %     mwPhase        double
% %                    detection phase (0 = absorption, pi/2 = dispersion)
% % 
% %     Temperature    double
% %                    temperature, in K
%
%
%   Opt: simulation options
%     chkCon         if equal to 1, after the first nSteps of the 
%                    trajectories are calculated, both inter- and intra-
%                    trajectory convergence is checked using the Gelman-
%                    Rubin R statistic such that R<1.1, and if this 
%                    condition is not satisfied, then propagation will be 
%                    extended by either a length of time equal to the 
%                    average of tcorr or by 20% more time steps, whichever 
%                    is larger
%
%     specCon        if equal to 1, after the first nOrients of the FID
%                    are calculated, both inter- and intra-FID convergence 
%                    are checked using the Gelman-Rubin R statistic such 
%                    that R<1.1, and if this condition is not satisfied, 
%                    then nOrients will be increased by 20% to simulate
%                    additional FIDs until R<1.1 is achieved
%
%     Verbosity      0: no display, 1: show info
%
%     Method         string
%                    Nitroxide: propagate the density matrix using an 
%                      analytical expression for the matrix exponential in 
%                      the m_S=-1/2 subspace (14N nitroxides only, faster)
%                    ISTOs: propagate the density matrix using
%                      irreducible spherical tensor operators (general, slower)
%
%    FFTWindow       1: use a Hamming window (default), 0: no window
%
%    truncate        double
%                    time point (in nanoseconds) at which to stop using 
%                    full quantum dynamics propagator and begin using an
%                    approximate propagator using correlation functions
%
%
%
%   MD: structure with molecular dynamics simulation parameters
%
%     dt             double
%                    time step (in s) for saving MD trajectory snapshots
%
%     tScale         double (optional)
%                    scale the time step of the simulation based on
%                    incorrect diffusion coefficients, e.g. molecules 
%                    solvated in TIP3P water have diffusion coefficients
%                    that are ~2.5x too large
%
%     tLag           double
%                    time lag (in s) for sampling the MD trajectory to 
%                    determine states and transitions, used for a Markov
%                    state model
%
%     DiffGlobal     double (optional)
%                    Diffusion coefficient for isotropic global rotational
%                    diffusion (s^-1)
%
%     TrajUsage      string (optional)
%                    Explicit: (default) use molecular orientations in
%                      trajectories directly as input for simulating the
%                      spectrum
%                    Resampling: coarse grain the trajectories by using the
%                      Euler angle probability distribution (for 
%                      pseudopotential) and orientational correlation 
%                      functions (for diffusion tensor) to perform further 
%                      stochastic rotational dynamics simulations
%
%     removeGlobal   integer
%                    1: (default) remove protein global diffusion
%                    0: no removal (e.g. if protein is fixed)
% 
%
%   Output:
%     B              numeric, size = (2*nSteps,1) 
%                    magnetic field (mT)
%
%     spc            numeric, size = (2*nSteps,1)
%                    derivative EPR spectrum
%
%     ExpectVal      numeric, size = (2*nSteps,1)
%                    expectation value of complex magnetization, 
%                    \langle S_{+} \rangle
%
%     t              numeric, size = (2*nSteps,1)
%                    simulation time axis (in s)

%    Opt.debug       struct
%                    various debugging options used for testing by
%                    developers
%
%                    EqProp: for computing the equilibrium propagator, 
%                    set to "time" (default) to only average over the time 
%                    axis of the Hamiltonian, yielding an approximate 
%                    propagator for each trajectory; set to "all" to 
%                    average over both time and trajectory axes, yielding a
%                    single approximate propagator to be used on the
%                    trajectory-averaged density matrix

%                    Markov: coarse grain the trajectories by using the
%                      side chain dihedral angles to form a Markov state 
%                      model


function varargout = cardamom(Sys,Exp,Par,Opt,MD)

% Preprocessing
% -------------------------------------------------------------------------

switch nargin
  case 0
    help(mfilename); return;
  case 2 % Par, Opt, and MD not specified, initialize them
    Par = [];
    Opt = [];
    MD = [];
  case 3 % Opt and MD not specified, initialize them
    Opt = [];
    MD = [];
  case 4 % MD not specified
    MD = [];
  case 5 % Sys, Par, Exp, Opt, MD specified
  otherwise
    error('Incorrect number of input arguments.')
end

error(chkmlver);

switch nargout
  case 0 % plotting
  case 1 % spc
  case 2 % B,spc
  case 3 % B,spc,ExpectVal
  case 4 % B,spc,ExpectVal,t
  otherwise
    error('Incorrect number of output arguments.');
end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end

% --------License ------------------------------------------------
LicErr = 'Could not determine license.';
Link = 'epr@eth'; eschecker; error(LicErr); clear Link LicErr
% --------License ------------------------------------------------

global EasySpinLogLevel;
global reverseStr
EasySpinLogLevel = Opt.Verbosity;


% Check Sys
% -------------------------------------------------------------------------

if ~isfield(Sys,'Nucs')
  Sys.Nucs = '14N';
  logmsg(0,'-- List of nuclei Sys.Nucs not specified. Using Sys.Nucs=''14N'' for a nitroxide spin label.');
end

[Sys,err] = validatespinsys(Sys);
error(err);

if Sys.nElectrons>1, error('cardamom does not support more than one electron.'); end
if Sys.nNuclei>1, error('cardamom does not support more than one nucleus.'); end

if isfield(Sys, 'lw')
  isBroadening = any(Sys.lw>0);
else
  isBroadening = 0;
end

% Check MD
% -------------------------------------------------------------------------

useMD = ~isempty(MD);

if useMD
  
  % scale the time axis if desired
  if isfield(MD,'tScale')
    tScale = MD.tScale;
  else
    tScale = 1;
  end
  
  if ~isfield(MD,'dt')
    error('The MD trajectory time step MD.dt must be given.')
  end
  MD.dt = tScale*MD.dt;
  
  if ~isfield(MD,'DiffGlobal')
    MD.DiffGlobal = [];
  end
  
  if ~isfield(MD,'nSteps')
    error('The number of MD trajectory time steps MD.nSteps must be given.')
  end
  
  % check type of MD trajectory usage
  if isfield(MD,'TrajUsage')
    if ~ischar(MD.TrajUsage)
      error('MD.TrajUsage must be a string.')
    end
    if ~any(strcmp({MD.TrajUsage},{'Explicit','Resampling'}))%,'Markov'}))
      errmsg = sprintf('Entry ''%s'' for MD.TrajUsage not recognized.', MD.TrajUsage);
      error(errmsg)
    end
  else
    MD.TrajUsage = 'Explicit';
  end
    
  if ~isfield(MD,'FrameTraj')
    error('The spin label frame trajectory MD.FrameTraj must be given.')
  end

  if ~isequal(size(MD.FrameTraj),[3,3,1,MD.nSteps])
    error('Frame trajectory must be of size (3,3,1,MD.nSteps).')
  end
  
  if ~isfield(MD,'removeGlobal')
    MD.removeGlobal = 1;
  end
  
  if MD.removeGlobal
    MD.RTraj = MD.FrameTrajwrtProt;
  else
    MD.RTraj = MD.FrameTraj;
  end

  % Check for orthogonality of rotation matrices
  RTrajInv = permute(MD.RTraj,[2,1,3,4]);

  if ~allclose(multimatmult(MD.RTraj,RTrajInv),...
               repmat(eye(3),1,1,size(MD.RTraj,3),size(MD.RTraj,4)),...
               1e-13)
    error('Rotation matrices obtained from frame trajectory are not orthogonal.')
  end

  MD.nTraj = size(MD.RTraj,3);  % this is assumed to be one for now
  
%   if strcmp(MD.TrajUsage,'Markov')
%     logmsg(1,'-- building Markov state model -----------------------------------------');
% 
%     if ~isfield(MD,'tLag')
%       error('If using a Markov model, the sampling lag time MD.tLag must be given.')
%     end
% 
%     % use lag time to sample the trajectory of dihedral angles such that 
%     % that the result is Markovian
%     nLag = ceil(MD.tLag/MD.dt);
%     MD.dihedrals = MD.dihedrals(:,:,1:nLag:end);
%     MD.dihedrals = permute(MD.dihedrals,[3,1,2]);
% 
%     % set the Markov chain time step based on the sampling lag time
%     Par.dt = tScale*MD.tLag;
%     
%     % we will need only the states when calling stochtraj_jump later
%     Opt.statesOnly = 1;
%     
%     if ~isfield(MD,'nStates')
%       MD.nStates = 48;
%     end
% 
%     % Perform k-means clustering
%     [MD.stateTraj,centroids] = clusterDihedrals(MD.dihedrals,MD.nStates,Opt.Verbosity);
%     MD.nSteps = size(MD.stateTraj, 1);  % TODO: find a way to process different step sizes here
%     
%     % remove chi3, as its dynamics are very slow on the typical MD
%     % timescale
%     MD.dihedrals = MD.dihedrals(:,logical([1,1,0,1,1]));
% 
%     % initialize HMM using clustering results
%     mu0 = centroids.';
%     for iState = 1:MD.nStates
%       idxState = MD.stateTraj==iState;
%       Sigma0(:,:,iState) = cov(MD.dihedrals(idxState,:));
%     end
% %     prior0 = normalise(rand(MD.nStates,1));
% %     transmat0 = mk_stochastic(rand(MD.nStates,MD.nStates));
% 
%     mixmat0 = ones(MD.nStates,1);
% %     mixmat0 = [];
% 
%     randints = randi(size(MD.stateTraj,1), MD.nStates, 1);
%     prior0 = MD.stateTraj(randints);
%     transmat0 = calc_TPM(MD.stateTraj,MD.nStates).';
%     
% %     % check if transmat0 is positive definite
% %     [~,isNotPosDef] = chol(transmat0);
% %     if isNotPosDef
% %       transmat0 = sqrtm(transmat0'*transmat0);
% %     end
% 
%     checkEmptyTrans = 1;
%     ProbRatioThresh = 1e-3;  % threshold in probability ratio for finding 
%                              % rarely visited states
%     while checkEmptyTrans
%       % run expectation-maximization algorithm on HMM model parameters
% %       ModelIn.prior = prior0;
% %       ModelIn.transmat = transmat0;
% %       ModelIn.mu = mu0;
% %       ModelIn.Sigma = Sigma0;
% % 
% %       [~, ModelOut] = cardamom_emghmm(MD.dihedrals, ModelIn, 20);
% % 
% %       prior1 = ModelOut.prior;
% %       transmat1 = ModelOut.transmat;
% %       mu1 = ModelOut.mu;
% %       Sigma1 = ModelOut.Sigma;
%       
% [~,prior1,transmat1,mu1,Sigma1,mixmat1] = ...
%     cardamom_emghmm(MD.dihedrals.',prior0,transmat0,mu0,Sigma0,mixmat0,...
%                     'max_iter',20,'verbose',Opt.Verbosity);
% 
%       nStates = size(transmat1,1);
%       EmptyTransList = zeros(nStates,1);
%       for iState = 1:nStates
%         % check for rare states
%         maxProb = max(transmat1(:));
%         StateTransProbs = [transmat1(iState,:).'; transmat1(:,iState)];
%         EmptyTransList(iState) = all(StateTransProbs./maxProb < ProbRatioThresh);
%       end
%       if any(EmptyTransList)
%         % reset model vars by removing entries for rarely visited states
%         idxNonEmpty = ~EmptyTransList;
%         prior0 = prior1(idxNonEmpty);
%         transmat0 = transmat1(idxNonEmpty,:);
%         transmat0 = transmat0(:,idxNonEmpty);
%         mu0 = mu1(:,idxNonEmpty);
%         Sigma0 = Sigma1(:,:,idxNonEmpty);
%         mixmat0 = mixmat1(idxNonEmpty);
%       else
%         % no rare states found, reset nStates if it has changed
%         checkEmptyTrans = 0;
%         MD.nStates = size(transmat1,1);
%       end
%     end
% %     prior1 = prior0;
% %     transmat1 = transmat0;
%   end
  % estimate rotational diffusion time scale
  FrameAcorr = autocorrfft(squeeze(MD.FrameTraj.^2), 2, 2, 1);

  N = round(MD.nSteps/4);

  % calculate correlation time
  time = linspace(0, N*MD.dt, N);
  tauR = max(cumtrapz(time,FrameAcorr(1:N),2),[],2);
%   [k,c,yfit] = exponfit(time, acorr(1:N), 2, 'noconst');
%   tauR = 1/max(k);

  tauR = mean(tauR);
  DiffLocal = 1/6/tauR;
  MD.tauR = tauR;

  clear RTrajInv

end

% Check local dynamics models
if ~isfield(Par,'Model')
  if useMD
    Par.Model = 'MD';
  else
    Par.Model = 'stochastic';
  end
end
if ~isempty(Par.Model)
  if ~strcmp(Par.Model,'stochastic') && ~strcmp(Par.Model,'jump') && ~strcmp(Par.Model,'MD')
    error('Model ''%s'' in Par.Model not recognized.',Par.Model);
  end
end

% Check Exp
% -------------------------------------------------------------------------

[Exp,FieldSweep,CenterField,CenterFreq,Sweep] = validate_exp('cardamom',Sys,Exp);

% if ~FieldSweep
%   error('cardamom does not support frequency sweeps.');  % TODO expand usage to include frequency sweep
% end

if FieldSweep
  omega0 = 2*pi*Exp.mwFreq*1e9;  % GHz -> rad s^-1
else
  omega0 = 2*pi*CenterFreq*1e9;  % MHz -> rad s^-1
end

if isfield(Exp,'SampleOrientation')
  SampleOrientation = Exp.SampleOrientation;
else
  SampleOrientation = [];
end


% set up horizontal sweep axis
if FieldSweep
  xAxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);  % field axis, mT
else
  xAxis = linspace(Exp.mwRange(1),Exp.mwRange(2),Exp.nPoints);  % field axis, GHz
end

% Check dynamics and ordering
% -------------------------------------------------------------------------

if useMD
  isDiffSim = strcmp(MD.TrajUsage,'Resampling');
elseif strcmp(Par.Model,'jump')
  isDiffSim = false;
else
  isDiffSim = true;
end

% FieldSweep = true;
Dynamics = validate_dynord('cardamom',Sys,FieldSweep,isDiffSim);

if isDiffSim
  DiffLocal = Dynamics.Diff;
  Sys.Diff = DiffLocal;
  tcorr = 1./6./DiffLocal;
end

if useMD
  Dynamics.DiffGlobal = MD.DiffGlobal;
end

if isfield(Sys,'Potential')
  isLocalPotential = true;
  LocalPotential = Sys.Potential;
else
  isLocalPotential = false;
end

% Check Par
% -------------------------------------------------------------------------

% Supply defaults
% default number of trajectories
if ~isfield(Par,'nTraj')&&useMD==0
  Par.nTraj = 100; 
  logmsg(0,'-- Number of trajectories Par.nTraj not given. Using %d trajectories.', Par.nTraj);
end
% Supply defaults
if ~isfield(Par,'nTraj')&&useMD==0, Par.nTraj = 100; end

if isfield(Par,'t')
  % time axis is given explicitly
  t = Par.t;
  tMax = max(t);
  nStepsQuant = numel(t);
  Par.Dt = t(2) - t(1);
  Par.dt = Par.Dt;
  nStepsStoch = round(tMax/Par.dt);
  % check for linearly spaced time axis
  absdev = abs(t/Par.Dt-(0:nStepsQuant-1));
  if max(absdev)>1e-13
    error('t does not appear to be linearly spaced.');
  end
  
elseif isfield(Par,'nSteps') && isfield(Par,'dt')
    % number of steps and time step are given
    nStepsQuant = Par.nSteps;
    if ~isfield(Par,'Dt')
      Par.Dt = Par.dt;
    end
    if Par.Dt<Par.dt
      error('The stochastic time step Par.dt must be less than or equal to Par.Dt.')
    end
    tMax = nStepsQuant*Par.Dt;
    nStepsStoch = round(tMax/Par.dt);

elseif isfield(Par,'nSteps') && isfield(Par,'tMax')
  % number of steps and max time are given
  tMax = Par.tMax;
  nStepsQuant = Par.nSteps;
  Par.Dt = tMax/Par.nSteps;
  Par.dt = Par.Dt;
  nStepsStoch = round(tMax/Par.dt);

else
  if isDiffSim
    Par.dt = min(tcorr)/10;
  elseif strcmp(Par.Model,'jump') && ~isfield(Par,'dt')
    error('The time step Par.dt must be specified when using an jump model.')
  elseif strcmp(Par.Model,'MD') && ~isfield(Par,'dt')
    error('The time step Par.dt must be specified when using an MD model.')
  end
  Par.Dt = Par.dt;
  logmsg(0,'-- No time parameters given. Using time step of %0.5g s.', Par.dt);
  if isfield(Par,'nSteps')
    nSteps = Par.nSteps;
  else
    nSteps = ceil(250e-9/Par.dt);
    logmsg(0,'-- Number of time steps not given. Using %d steps.', nSteps);
  end
  nStepsStoch = nSteps;
  nStepsQuant = nSteps;
end

dtQuant = Par.Dt;
dtStoch = Par.dt;

% decide on a simulation model based on user input
if useMD
  if ~isfield(Par,'Model')
    % no Model given
    Par.Model = 'MD';
  elseif ~strcmp(Par.Model,'MD')
    error('Mixing stochastic simulations with MD simulation input is not supported.')
  end
  
  if MD.nTraj > 1, error('Using multiple MD trajectories is not supported.'); end
  
  % determine if time block averaging is to be used
  if strcmp(MD.TrajUsage,'Explicit')
    % check MD.dt
    if Par.Dt<MD.dt
      error('Par.Dt must be greater than MD.dt.')
    end
    Par.isBlock = 1;
  else
    % check Par.Dt
    if Par.dt<Par.Dt
      Par.isBlock = 1;
    else
      % same step size
      Par.isBlock = 0;
    end
  end
  
  % find block properties for block averaging
  if Par.isBlock
    if strcmp(MD.TrajUsage,'Explicit')
      [Par.nBlocks,Par.BlockLength] = findblocks(Par.Dt, MD.dt, MD.nSteps);
    else
      [Par.nBlocks,Par.BlockLength] = findblocks(Par.Dt, Par.dt, nStepsStoch);
    end
  end
  
  % process single long trajectory into multiple short trajectories
  if strcmp(MD.TrajUsage,'Explicit')
%     Par.lag = ceil(3*MD.tauR/Par.Dt);  % use 2 ns lag between windows
    Par.lag = ceil(2e-9/Par.Dt);  % use 2 ns lag between windows
    if Par.nSteps<Par.nBlocks
      % Par.nSteps not changed from user input
      Par.nTraj = floor((Par.nBlocks-Par.nSteps)/Par.lag) + 1;
    else
      Par.nSteps = Par.nBlocks;
      Par.nTraj = 1;
    end
  end
  
else
  % no MD simulation trajectories provided, so perform stochastic dynamics
  % simulations internally
  
  if ~isfield(Par,'Model')
    Par.Model = 'stochastic';
  end
  
  % check for time coarse-graining (time block averaging)
  if Par.dt<Par.Dt
    Par.isBlock = 1;
    [Par.nBlocks,Par.BlockLength] = findblocks(Par.Dt, Par.dt, nStepsStoch);
  else 
    % same step size
    Par.isBlock = 0;
  end
    
end

LocalDynamicsModel = Par.Model;


% Check Opt
% -------------------------------------------------------------------------

if ~isfield(Opt,'Method')
  Opt.Method = 'Nitroxide';
end

if ~isfield(Opt,'truncate')
  Opt.truncate = 0;
end

if Opt.truncate && strcmp(Opt.Method,'Nitroxide')
  error('Correlation function propagation is only available for ISTOs method.')
end

if isfield(Opt,'FFTWindow')
  FFTWindow = Opt.FFTWindow;
else
  FFTWindow = 1;
end

if isfield(Opt,'specCon')
  specCon = Opt.specCon;
else
  specCon = 0;
end

if ~isfield(Opt,'debug')
  Opt.debug = [];
end

% Debugging options
if ~isfield(Opt,'debug')
  Opt.debug.EqProp = 'all';
else
  if ~isfield(Opt.debug,'EqProp')
    Opt.debug.EqProp = 'all';
  else
    switch Opt.debug.EqProp
      case 'time'
        Opt.debug.EqProp = 'time';
      case 'all'
        Opt.debug.EqProp = 'all';
      otherwise
        error('Opt.debug.EqProp value not recognized.')
    end
  end
end


% Check model for local diffusion
% -------------------------------------------------------------------------
switch LocalDynamicsModel
  case 'stochastic'
    
    DiffLocal = Dynamics.Diff;
  
  case 'jump'
    
  case 'MD' % TODO process RTraj based on size of input
    
    if ~isfield(Par,'nOrients')
      error('nOrients must be specified for the MD model.')
    end
    
%     if strcmp(Opt.Method, 'ISTOs')
%       % this method uses quaternions, not rotation matrices, so convert
%       % MD.RTraj to quaternions here before the simulation loop
    RTrajLocal = MD.RTraj;
    qTrajLocal = rotmat2quat(MD.RTraj);
    
    if ~strcmp(MD.TrajUsage,'Explicit')
      switch MD.TrajUsage
        case 'Resampling'
          M = size(MD.FrameTraj,4);

          % calculate orienting potential energy function
          theta = squeeze(acos(MD.FrameTraj(3,3,:,1:M)));
          phi = squeeze(atan2(MD.FrameTraj(3,2,:,1:M), MD.FrameTraj(3,1,:,1:M)));
          psi = squeeze(atan2(-MD.FrameTraj(2,3,:,1:M), MD.FrameTraj(1,3,:,1:M)));
%           theta = squeeze(acos(MD.FrameZ(3,:,:,1:M)));
%           phi = squeeze(atan2(MD.FrameY(3,:,:,1:M), MD.FrameX(3,:,:,1:M)));
%           psi = squeeze(atan2(-MD.FrameZ(2,:,:,1:M), MD.FrameZ(1,:,:,1:M)));

          phi = phi + 2*pi*(phi<0);
          psi = psi + 2*pi*(psi<0);

          nBins = 90;
          phiBins = linspace(0, 2*pi, nBins);
          thetaBins = linspace(0, pi, nBins/2);
          psiBins = linspace(0, 2*pi, nBins);

          [pdf, ~] = histcnd([phi,theta,psi], {phiBins,thetaBins,psiBins});

          pdf(end,:,:) = pdf(1,:,:);  % FIXME why does it truncate to zero in the phi direction?
          pdf = smooth3(pdf,'gaussian');
          pdf(pdf<1e-14) = 1e-14;  % put a finite floor on histogram
          Sys.Potential = -log(pdf);
%           pdf = smoothn(pdf);
%         case 'Markov'
%           
%           RTrajLocal = RTrajLocal(:,:,1,1:nLag:end);
%           qTrajLocal = rotmat2quat(RTrajLocal);
          
      end
      
    end
    
  % these variables could be huge and are no longer needed, so delete them 
  % now
  MD = rmfield(MD, 'FrameTraj');
  MD = rmfield(MD, 'RTraj');
%   MD = rmfield(MD, 'FrameTrajwrtProt');
%   MD = rmfield(MD, 'RProtDiff');
    
end

% Generate grids for powder averaging in the lab frame
% -------------------------------------------------------------------------

% default number of orientations
if isfield(Par,'nOrients')
  nOrients = Par.nOrients;
else
  nOrients = Par.nTraj;
  logmsg(0,'-- Par.nOrients not specified. Using %d orientations.', nOrients);
end

Orients = [];
if isfield(Par,'Orients')
  Orients = Par.Orients;
  if isvector(Orients)
    % assure that Orients is a column vector
    Orients = Orients(:);
    if nOrients>1
      % if only one orientation, then repeat it nOrients times
      Orients = repmat(Orients,[1,nOrients]);
    end
  end
end

% assign lab orientations for powder averaging
if ~isempty(Orients)
  gridPhi = Orients(1,:);
  gridTheta = Orients(2,:);
else
  if specCon
    % generate Sobol sequence over a spiral grid
    nOrients = ceil(nOrients/2);
    skip = 0;
    gridPts = 2*cardamom_sobol_generate(1,nOrients,skip)-1;
    gridPhi = sqrt(pi*nOrients)*asin(gridPts);
    gridTheta = acos(gridPts);
  else
    % generate a spiral grid
    gridPts = linspace(-1,1,nOrients);
    gridPhi = sqrt(pi*nOrients)*asin(gridPts);
    gridTheta = acos(gridPts);
  end
end

logmsg(1,'-- time domain simulation -----------------------------------------');
logmsg(1,'-- Model: %s -----------------------------------------', LocalDynamicsModel);
logmsg(1,'-- Method: %s -----------------------------------------', Opt.Method);

% Run simulation
% -------------------------------------------------------------------------

clear cardamom_propagatedm

converged = 0;
iter = 1;
spcArray = [];
nOrientsTot = 0;

t = linspace(0, nStepsQuant*Par.Dt, nStepsQuant);

while ~converged
  tic

  % trajectories might differ in length, so we need cells for allocation
  ExpectVal = {};
  tCell = [];
  reverseStr = [];

  % temporary cells to store intermediate results
  iExpectVal = cell(1,nOrients);
  itCell = cell(1,nOrients);

  for iOrient = 1:nOrients
    % Generate trajectory of local dynamics
    switch LocalDynamicsModel
      
      case 'stochastic'
        
        Sys.Diff = Dynamics.Diff;
        Par.dt = dtStoch;
        Par.nSteps = nStepsStoch;
        if isLocalPotential
          Sys.Potential = LocalPotential;
          Par.nSteps = 2*nStepsStoch;
        end
        [~, RTrajLocal, qTrajLocal] = stochtraj_diffusion(Sys,Par,Opt);
        if isLocalPotential
          RTrajLocal = RTrajLocal(:,:,:,nStepsStoch+1:end);
          qTrajLocal = qTrajLocal(:,:,nStepsStoch+1:end);
        end
        
      case 'jump'
        
        Par.dt = dtStoch;
        Par.nSteps = nStepsStoch;
        [~, RTrajLocal, qTrajLocal] = stochtraj_jump(Sys,Par,Opt);
        
      case 'MD'
        
        switch MD.TrajUsage
          case 'Explicit'
            
            % the MD trajectories are not changing, so RTraj and qTraj were
            % processed earlier outside of the loop
            
          case 'Resampling'
            
            Sys.ProbDensFun = pdf;
            Sys.Diff = DiffLocal;
            if isLocalPotential
              Sys.Potential = LocalPotential;
            end
            Par.dt = dtStoch;
            Par.nSteps = nStepsStoch;
            [~, RTrajLocal, qTrajLocal] = stochtraj_diffusion(Sys,Par,Opt);
            
%           case 'Markov'
%             
%             Sys.States0 = rejectionsample(MD.nStates, prior1, Par.nTraj);
%             Sys.TransProb = transmat1.';
%             Par.dt = dtStoch;
%             Par.nSteps = nStepsStoch;
%             Opt.statesOnly = true;
%             [~, stateTraj] = stochtraj_jump(Sys,Par,Opt);
%             
%             Par.stateTraj = stateTraj;
            
        end
        
    end
    
    if strcmp(Opt.Method,'ISTOs')
      Par.qTraj = qTrajLocal;
    else
      Par.RTraj = RTrajLocal;
    end
    
    % Generate trajectory of global dynamics
    includeGlobalDynamics = ~isempty(Dynamics.DiffGlobal);
    if includeGlobalDynamics
      Sys.Diff = Dynamics.DiffGlobal;
     if isLocalPotential
       Sys.Potential = LocalPotential;
     end
      Par.dt = dtQuant;
      Par.nSteps = nStepsQuant;
      [~, ~, qTrajGlobal] = stochtraj_diffusion(Sys,Par,Opt);
    end
    
    % Combine global trajectories with starting orientations
    qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
                  [1,Par.nTraj,nStepsQuant]);
    if includeGlobalDynamics
      qLab = quatmult(qLab,qTrajGlobal);
    end
    
    if strcmp(Opt.Method,'ISTOs')
      Par.qLab = qLab;
    else
      Par.RLab = quat2rotmat(qLab);
    end
    
    % propagate the density matrix
    Par.nSteps = nStepsQuant;
    Par.Dt = dtQuant;
    Par.dt = dtStoch;
    Sprho = cardamom_propagatedm(Sys,Par,Opt,MD,omega0,CenterField);
    
    % average over trajectories
    if strcmp(Opt.debug.EqProp,'time')
      Sprho = squeeze(mean(Sprho,3));
    end
    
    iExpectVal{1,iOrient} = 0;
    % calculate the expectation value of S_{+}
    for k = 1:size(Sprho,1)
      iExpectVal{1,iOrient} = iExpectVal{1,iOrient} + squeeze(Sprho(k,k,:));  % take traces TODO try to speed this up using equality tr(A*B)=sum(sum(A.*B))
    end
    
    if Opt.Verbosity
      updateuser(iOrient,nOrients)
    end
    
    itCell{1,iOrient} = t.';
    
    iOrient = iOrient + 1;
    
  end

  % Store simulations at new starting orientations from each iteration
  ExpectVal = cat(2, ExpectVal, iExpectVal);
  tCell = cat(2, tCell, itCell);

  % Perform FFT
  % -------------------------------------------------------------------------

  % windowing
  if FFTWindow
    hamming = cellfun(@(x) 0.54 + 0.46*cos(pi*x/max(x)), tCell, 'UniformOutput', false);
    hann = cellfun(@(x) 0.5*(1 + cos(pi*x/max(x))), tCell, 'UniformOutput', false);
    Window = hann;
    ExpectVal = cellfun(@times, ExpectVal, Window, 'UniformOutput', false);
  end

  % zero padding for FFT to ensure sufficient B-field resolution (at most 0.1 G)
  % expectval = cell2mat(cellfun(@(x) zeropad(x, maxlength), expectval, 'UniformOutput', false));
  % tlong = (0:Par.dt:maxlength*Par.dt);

  Bres = 0.1; % G
  tReq = 1/(mt2mhz(Bres/10)*1e6); % mT -> s

  % if max(t)<treq
  tMax = max(cellfun(@(x) max(x), tCell));
  if tMax<tReq
    M = ceil(tReq/Par.Dt);
  else
    M = ceil(tMax/Par.Dt);
  end
  ExpectVal = cell2mat(cellfun(@(x) zeropad(x, M), ExpectVal, 'UniformOutput', false));
  tLong = linspace(0, M*Par.Dt, M).';

  % convolute for linewidth broadening
  if isBroadening
    if Sys.lw(1)>0
      % Gaussian broadening
      w = mt2mhz(Sys.lw(1))*1e6;  % FWHM in Hz
      alpha = pi^2*w^2/(4*log(2));
      ExpectVal = bsxfun(@times,exp(-alpha*tLong.^2),ExpectVal);
    end
    if numel(Sys.lw)==2
      % Lorentzian broadening
      TL = Dynamics.T2; 
      ExpectVal = bsxfun(@times,exp(-tLong/TL),ExpectVal);
    end
  end
  
  % Multiply by t for differentiation and take the FFT
  spcArray = cat(2, spcArray, imag(fftshift(fft(bsxfun(@times,ExpectVal,tLong), [], 1))));
  spcNew = mean(spcArray,2);

  if specCon
    if iter==1
      spcLast = spcNew;
      rmsdNew = NaN;
    else
      span = max(spcNew)-min(spcNew);
      rmsdNew = sqrt(mean((spcNew-spcLast).^2))/span;
      if rmsdNew<2e-2
        converged = 1;
      elseif iter>2
        rmsdPctChange = abs(100*(rmsdNew-rmsdLast)/rmsdLast)
        converged = rmsdPctChange<50;
      end
    end

  else
    converged = 1;
  end
  
  if ~converged
    % double the number of orientations in powder averaging
    msg = sprintf('Convergence not achieved. Propagation is being extended.\n');
    if Opt.Verbosity
      fprintf(msg);
    end
    rmsdLast = rmsdNew;
    spcLast = spcNew;  % store for comparison after next iteration completes
    nOrientsTot = nOrientsTot + nOrients;
    skip = iter*nOrientsTot;  % seed Sobol sequence generator for next iteration
    
    nOrients = nOrientsTot;
    if ~isempty(Orients)
      gridPhi = repmat(Orients(1,:),[1,2^(iter-1)]);
      gridTheta = repmat(Orients(2,:),[1,2^(iter-1)]);
    else
      gridPts = 2*cardamom_sobol_generate(1,nOrients,skip)-1;
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    end
    iter = iter + 1;
  else
    spcAvg = spcNew;
    minsTot = floor(toc/60);
    if Opt.Verbosity
      msg = sprintf('Done!\nTotal simulation time: %d:%2.0f\n',minsTot,mod(toc,60));
      fprintf(msg);
    end
  end

end

clear RTraj RTrajInv qTraj
% Par = rmfield(Par,{'RTraj','qTraj'});

if strcmp(LocalDynamicsModel, 'Molecular Dynamics')
  % these variables can take up a lot of memory and might prevent the user 
  % from implementing a fine enough grid for powder averaging 
  clear MD
end

freq = 1/(Par.Dt*M)*(-M/2:M/2-1);  % TODO check for consistency between FieldSweep and FFT window/resolution

% center the spectrum around the isotropic component of the g-tensor
if FieldSweep
  fftAxis = mhz2mt(freq/1e6+Exp.mwFreq*1e3, mean(Sys.g));  % MHz -> mT, note use of g0, not ge
else
  fftAxis = freq/1e9+mt2mhz(Exp.Field,mean(Sys.g))/1e3;  % GHz
end

% interpolate over horizontal sweep range
if FieldSweep
  outspc = interp1(fftAxis, spcAvg, xAxis);
else
  spcAvg = spcAvg(end:-1:1);   % reverse the axis for frequency sweep
  outspc = interp1(fftAxis, spcAvg, xAxis);
  outspc = cumtrapz(xAxis(end:-1:1),outspc);  % frequency sweeps outputs the absorption
end

% average over trajectories for expectation value output
ExpectVal = mean(ExpectVal, 2);
ExpectVal = ExpectVal(1:Par.nSteps);

% Final processing
% -------------------------------------------------------------------------

switch (nargout)
case 0
  cla
  if FieldSweep  % TODO fix output plotting
    if (xAxis(end)<10000)
      plot(xAxis,outspc);
      xlabel('magnetic field (mT)');
    else
      plot(xAxis/1e3,outspc);
      xlabel('magnetic field (T)');
    end
    axis tight
    ylabel('intensity (arb.u.)');
    title(sprintf('%0.8g GHz, %d points',Exp.mwFreq,numel(xAxis)));
  else
    if (xAxis(end)<1)
      plot(xAxis*1e3,spc);
      xlabel('frequency (MHz)');
    else
      plot(xAxis,spc);
      xlabel('frequency (GHz)');
    end
    axis tight
    ylabel('intensity (arb.u.)');
    title(sprintf('%0.8g mT, %d points',Exp.Field,numel(xAxis)));
  end
case 1
  varargout = {outspc};
case 2
  varargout = {xAxis,outspc};
case 3
  varargout = {xAxis,outspc,ExpectVal};
case 4
  varargout = {xAxis,outspc,ExpectVal,t};
end

clear global EasySpinLogLevel
    
end

% Helper functions
% -------------------------------------------------------------------------



function [stateTraj,centroids] = clusterDihedrals(dihedrals,nStates,verbosity)

chi1 = dihedrals(:,1);
chi2 = dihedrals(:,2);
% chi3 = dihedrals(:,3);
chi4 = dihedrals(:,4);
chi5 = dihedrals(:,5);

dihedrals = [wrapTo2Pi(chi1), wrapTo2Pi(chi2), wrapTo2Pi(chi4), chi5];

% initialize cluster centroids
% chi1Min = wrapTo2Pi([-60;65;180]/180*pi);
% chi2Min = wrapTo2Pi([75;180]/180*pi);
% chi4Min = wrapTo2Pi([75;8;-100]/180*pi);
% chi5Min = wrapTo2Pi([180;77]/180*pi);

[stateTraj,centroids] = cardamom_kmeans(dihedrals, nStates, 20, verbosity);

end

function TPM = calc_TPM(stateTraj, nStates)

Nij = zeros(nStates);

stateLast = stateTraj(1);
nSteps = length(stateTraj);

for iStep = 2:nSteps
  stateNew = stateTraj(iStep);
  Nij(stateLast,stateNew) = Nij(stateLast,stateNew) + 1;
  stateLast = stateNew;
end

Nij = (Nij + Nij.')/2;
TPM = Nij./sum(Nij,1);

end

function updateuser(iOrient,nOrient)
% Update user on progress

% persistent reverseStr
global reverseStr

% if isempty(reverseStr), reverseStr = []; end

avgTime = toc/iOrient;
secsLeft = (nOrient - iOrient)*avgTime;
minsLeft = floor(secsLeft/60);

secsElap = toc;
minsElap =  floor(secsElap/60);

msg1 = sprintf('Lab orientation: %d/%d\n', iOrient, nOrient);
if avgTime<1.0
  msg2 = sprintf('%2.1f orientations/s\n', 1/avgTime);
else
  msg2 = sprintf('%2.1f s/orientation\n', avgTime);
end
msg3 = sprintf('Time elapsed: %d:%2.0f\n', minsElap, mod(secsElap,60));
msg4 = sprintf('Time remaining (predicted): %d:%2.0f\n', minsLeft, mod(secsLeft,60));
msg = [msg1, msg2, msg3, msg4];

fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));

end

function [nBlocks,BlockLength] = findblocks(Dt, dt, nSteps)
% time step of molecular simulation, dt, is usually much smaller than
% the quantum propagation time step, Dt, so determine the size of the 
% averaging block here (Dt/dt), then set nTraj and nSteps accordingly to 
% number of blocks and the size of the block, respectively

% size of averaging blocks
BlockLength = ceil(Dt/dt);

% size of trajectory after averaging, i.e. coarse-grained trajectory length
nBlocks = floor(nSteps/BlockLength);

BlockLength = BlockLength;
nBlocks = nBlocks;

end

function xSamples = rejectionsample(maxX, p, nSamples)

xGrid = 1:maxX;

c = max(p(:));

iSample = 1;
xSamples = zeros(1,nSamples);
while iSample < nSamples + 1
  xProposal = randi(maxX);
  q = c;
  if rand() < p(xProposal)/q
    xSamples(1,iSample) = xGrid(xProposal);
    iSample = iSample + 1;
  end
end

end

function  y = zeropad(x, M)
  N = length(x);
  if iscolumn(x), y = [x; zeros(M-N, 1)]; end
  if isrow(x), y = [x, zeros(1, M-N)]; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           %%%%%%%%%%%%
%           
%           bins = size(pdf, 1);
% 
%           abins = bins;  % use less points for easier heat map visualization
%           bbins = bins/2;
%           gbins = bins;
%           
%           alphaBins = linspace(-pi, pi, abins);
%           betaBins = pi-acos(linspace(-1, 1, bbins));
%           gammaBins = linspace(-pi, pi, gbins);
%           
%           for iTraj=1:Par.nTraj
%             % use a "burn-in method" by taking last half of each trajectory
%             [alpha, beta, gamma] = quat2euler(qTraj(:,iTraj,:),'active');
%           %   [alpha, beta, gamma] = quat2euler(qTraj(:,iTraj,N:end));
%             alpha = squeeze(alpha);
%             beta = squeeze(beta);
%             gamma = squeeze(gamma);
% 
%           %   q = squeeze(qTraj(:,iTraj,N:end));
%           %   [alpha, beta, gamma] = quat2angle(permute(q,[2,1]), 'ZYZ');
% 
%             % calculate 3D histogram using function obtained from Mathworks File Exchange
%             [Hist3D(:,:,:,iTraj),~] = histcnd([alpha,beta,gamma],...
%                                               {alphaBins,betaBins,gammaBins});
%           end
%           
%           hist_tauD100 = mean(Hist3D,4);
%           save('hist_tauD100.mat', 'hist_tauD100')
%           
%           %%%%%%%%%%%%%

% used for plotting PseudoPotFun on a sphere

% %       pad = (PseudoPotFun(1,:,:) + PseudoPotFun(end,:,:))/2;
%       
%       pdf = smooth3(pdf, 'gaussian');
% 
%       yy = permute(mean(pdf, 3), [2, 1]);
%       yy = yy/max(yy(:));
%       yy(:,end) = yy(:,1);
%       
%       theta = linspace(0, pi, size(yy,1));                   % polar angle
%       phi = linspace(0, 2*pi, size(yy,2));                   % azimuth angle
%       
%       [Phi, Theta] = meshgrid(phi, theta);
%       r = 1.0;
%       amplitude = 1.0;
% %       rho = r + amplitude*yy;
%       rho = yy;
%       
%       x = r.*sin(Theta).*cos(Phi);
%       y = r.*sin(Theta).*sin(Phi);
%       z = r.*cos(Theta);
% 
%       surf(x, y, z, rho, ...
%            'edgecolor', 'none', ...
%            'facecolor', 'interp');
%       % title('$\ell=0, m=0$')
% 
% 
% %       shading interp
% 
%       axis equal off      % set axis equal and remove axis
%       view(90,30)         % set viewpoint
%       set(gca,'CameraViewAngle',6);
%       
%       left = 0.5;
%       bottom = 0.5;
%       width = 4;     % Width in inches
%       height = 4;    % Height in inches
% 
%       set(gcf,'PaperUnits','inches');
%       set(gcf,'PaperSize', [8 8]);
%       set(gcf,'PaperPosition',[left bottom width height]);
%       set(gcf,'PaperPositionMode','Manual');
%       
%       print('ProbabilityDist', '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%           nBins = 50;
%           
%           Hist3D = zeros(nBins, nBins, nBins, Par.nTraj);
% 
%           % alphaBins = pi-2*acos(linspace(-1, 1, abins));
%           alphaBins = linspace(-pi, pi, nBins);
%           betaBins = pi-acos(linspace(-1, 1, nBins));
% %           betaBins = linspace(0, pi, nBins);
%           gammaBins = linspace(-pi, pi, nBins);
%           
% %           M = size(qTraj,3)/2;
%           M = 1;
%           
% 
%           for iTraj=1:Par.nTraj
%             % use a "burn-in method" by taking last half of each trajectory
% %             [alpha, beta, gamma] = quat2euler(qTraj(:,iTraj,M:end));
%             [alpha, beta, gamma] = quat2euler(qTraj(:,iTraj,M:end), 'active');
%             alpha = squeeze(alpha);
%             beta = squeeze(beta);
%             gamma = squeeze(gamma);
% 
% 
%             % calculate 3D histogram using function obtained from Mathworks File Exchange
%             [Hist3D(:,:,:,iTraj),~] = histcnd([alpha,beta,gamma],...
%                                               {alphaBins,betaBins,gammaBins});
%           end
% 
%           Hist3D = mean(Hist3D, 4);  % average over all trajectories
%           Hist3D = Hist3D/trapz(Hist3D(:));  % normalize
%           
%           save('Hist3D.mat', 'Hist3D')
%           
%           M = 1;
          
%           figure(1)
%           subplot(1,2,1)
%           slice(alphaBins, ...
%                 betaBins, ...
%                 gammaBins, ...
%                 permute(pdfDec, pidx), ...
%                 0, pi/2, 0)
%           xlabel('alpha')
%           ylabel('beta')
%           zlabel('gamma')
%           colormap hsv
%           
%           subplot(1,2,2)
%           slice(alphaBins, ...
%                 betaBins, ...
%                 gammaBins, ...
%                 permute(Hist3D, pidx), ...
%                 0, pi/2, 0)
%           xlabel('alpha')
%           ylabel('beta')
%           zlabel('gamma')
%           colormap hsv

  %         N = 50;
  %         Aedges = linspace(-pi/2, pi/2, N);
  %         Bedges = pi-acos(linspace(-1, 1, N));
  % 
  %         for iTraj=1:Par.nTraj
  %             [temp,~] = histcounts2(alpha(iTraj,:), beta(iTraj,:), Aedges, Bedges);
  %             Hist2D(:,:,iTraj) = temp;
  %             [temp2,~] = histcounts2(alpha(iTraj,:), beta(iTraj,:), Aedges, linspace(0,pi,N));
  %             Hist2D2(:,:,iTraj) = temp2;
  %           [temp,~] = histcnd([alpha(iTraj,:).',beta(iTraj,:).',gamma(iTraj,:).'],...
  %                                         {PhiBins.',ThetaBins.',PsiBins.'});
  %                                       
  %           Hist3D(:,:,:,iTraj) = permute(temp, [2, 1, 3]);
  %         end
  %         
  % % %       pad = (PseudoPotFun(1,:,:) + PseudoPotFun(end,:,:))/2;
  %         HistAvg = mean(Hist3D, 4);
  %         HistAvg(end,:,:) = HistAvg(1,:,:);
  % 
  % %         HistAvg = smooth3(HistAvg, 'gaussian');
  % 
  % %         yy = permute(mean(HistAvg, 3), [2, 1]);
  %         yy = mean(HistAvg, 3);
  %         yy = yy/max(yy(:));
  %         yy(:,end) = yy(:,1);
  % 
  %         theta = linspace(0, pi, size(yy,1));                   % polar angle
  %         phi = linspace(0, 2*pi, size(yy,2));                   % azimuth angle
  % 
  %         [Phi, Theta] = meshgrid(phi, theta);
  %         radius = 1.0;
  %         amplitude = 1.0;
  %   %       rho = radius + amplitude*yy;
  %         rho = yy;
  % 
  %         r = radius.*sin(Theta);    % convert to Cartesian coordinates
  %         x = r.*cos(Phi);
  %         y = r.*sin(Phi);
  %         z = radius.*cos(Theta);
  % 
  %         surf(x, y, z, rho, ...
  %              'edgecolor', 'none', ...
  %              'facecolor', 'interp');
  %         % title('$\ell=0, m=0$')
  % 
  % 
  %   %       shading interp
  % 
  %         axis equal off      % set axis equal and remove axis
  %         view(90,30)         % set viewpoint
  %         set(gca,'CameraViewAngle',6);
  %         
  %         HistTot = HistTot + mean(Hist3D,4);