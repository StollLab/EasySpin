% cardamom  Trajectory-based simulation of CW-EPR spectra.
%
%   cardamom(Sys,Exp,Par)
%   cardamom(Sys,Exp,Par,Opt)
%   cardamom(Sys,Exp,Par,Opt,MD)
%   spc = cardamom(...)
%   [B,spc] = cardamom(...)
%   [B,spc,TDSignal,t] = cardamom(...)
%   
%   Computes a CW-EPR spectrum of an S=1/2 spin label using stochastic or
%   molecular dynamics trajectories.
%   
%   Sys: stucture with spin system's static and dynamical parameters
%     tcorr          double or numeric vector, size = (1,3)
%                    correlation time (in seconds)
%     logtcorr       double or numeric vector, size = (1,3)
%                    log10 of rotational correlation time (in seconds)
%     Diff           double or numeric vector, size = (1,3)
%                    diffusion rate (rad^2 s^-1)
%     logDiff        double or numeric vector, size = (1,3)
%                    log10 of diffusion rate (rad^2 s^-1)
%
%         All fields can have 1 (isotropic), 2 (axial) or 3 (rhombic) elements.
%         Precedence: logtcorr > tcorr > logDiff > Diff.
%
%     DiffGlobal     double or numeric vector, size = (1,3)
%                    global diffusion rate (s^-1)
%     Potential      defines an orienting potential using one of the following:
%                    a) set of LMKs and lambdas for an expansion in terms of
%                       Wigner functions [L1 M1 K1 lambda1; L2 M2 K2 lambda2]
%                       lambda can be real or complex
%                    b) 3D array defining the potential over a grid, with the
%                       first dimension representing alpha [0,2*pi], the second
%                       beta [0,pi], and the third gamma [0,2*pi]
%                    c) a function handle for a function that takes three
%                       arguments (alpha, beta, and gamma) and returns the value
%                       of the orientational potential for that orientation. The
%                       function should be vectorized, i.e. work with arrays of
%                       alpha, beta, and gamma.
%     TransRates     numeric, size = (nStates,nStates)
%                    transition rate matrix describing inter-state dynamics
%                    for kinetic Monte Carlo simulations
%     TransProb      numeric, size = (nStates,nStates)
%                    transition probability matrix describing inter-state 
%                    dynamics for kinetic Monte Carlo simulations, note
%                    that a time step must be given to use Sys.TransProb
%                    (alternative input to TransRates; ignored if TransRates
%                    is given)
%     Orientations   numeric, size = (3,nStates)
%                    Euler angles for each state's orientation
%     lw             double or numeric vector, size = (1,2)
%                    vector with FWHM residual broadenings (in mT)
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
%     lwpp           peak-to-peak linewidths (mT), same format as lw
%                    use either lw or lwpp
%
%   Exp: experimental parameter settings
%     mwFreq         microwave frequency, in GHz (for field sweeps)
%     Range          sweep range, [sweepmin sweepmax], in mT (for field sweep)
%     CenterSweep    sweep range, [center sweep], in mT (for field sweep)
%     nPoints        number of points
%     Harmonic       detection harmonic: 0, 1 (default), 2
%     ModAmp         peak-to-peak modulation amplitude, in mT (field sweeps only)
%     mwPhase        detection phase (0 = absorption, pi/2 = dispersion)
%     Temperature    temperature, in K
%
%   Par: structure with simulation parameters
%     Model      model for spin label dynamics
%                'diffusion': Brownian rotation diffusion with given rotational
%                   diffusion tensor and ordering potential
%                'jump': Markovian jumps between a given set of discrete states
%                'MD-direct': use molecular orientations in MD
%                  trajectories directly as input for simulating the
%                  spectrum
%                'MD-HBD': coarse grain the MD trajectories by using 
%                  the Euler angle probability distribution (for 
%                  pseudopotential) from the spin label's orientations 
%                  to perform further stochastic rotational dynamics 
%                  simulations
%                'MD-HMM': coarse grain the MD trajectories by using 
%                  the spin label's side chain dihedral angles to build 
%                  a hidden Markov model model to perform further 
%                  stochastic jump dynamics simulations
%     dtSpatial  spatial dynamics propagation time step (in seconds)
%                (not used for 'MD-direct')
%     dtSpin     spin dynamics propagation time step (in seconds)
%     nSteps     number of time steps per simulation
%     nTraj      number of trajectories
%     OriStart   numeric, size = (3,1), (1,3), or (3,nTraj)
%                Euler angles for starting orientation(s)
%     nOrients   number of lab-to-molecule orientations to loop over
%     Orients    numeric matrix, size = (nOrients,2)
%                (optional) (phi,theta) angles of lab-to-molecule 
%                orientations. If not given, these are chosen as points
%                on a spherical spiral grid
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
%     specCon        if equal to 1, after the first nOrients of the FID
%                    are calculated, both inter- and intra-FID convergence 
%                    are checked using the Gelman-Rubin R statistic such 
%                    that R<1.1, and if this condition is not satisfied, 
%                    then nOrients will be increased by 20% to simulate
%                    additional FIDs until R<1.1 is achieved
%     Verbosity      0: no display, 1: show info
%     Method         string
%                    fast: propagate the density matrix using an 
%                      analytical expression for the matrix exponential in 
%                      the m_S=-1/2 subspace (S=1/2 with up to one nucleus)
%                    ISTOs: propagate the density matrix using
%                      irreducible spherical tensor operators (general, slower)
%     FFTWindow       1: use a Hamming window (default), 0: no window
%     nTrials        number of initialization trials for k-means
%                    clustering; used for the Markov method
%     LagTime        lag time for sliding window processing (only used for
%                    'MD-direct')
%
%   MD: structure with molecular dynamics simulation parameters
%
%     FrameTraj      frame trajectories, size (3,3,nSteps,nTraj)
%     FrameTrajwrtProtein frame trajecoties with respect to protein frame,
%                    size (3,3,nSteps,nTraj)
%     dt             time step (in s) for saving MD trajectory snapshots
%     tLag           time lag (in s) for sampling the MD trajectory to 
%                    determine states and transitions, used for the hidden
%                    Markov model
%     nStates        number of states in the hidden Markov model
%     DiffGlobal     diffusion coefficient for isotropic global rotational
%                    diffusion (s^-1)
%     removeGlobal   integer
%                    1: (default) remove protein global diffusion
%                    0: no removal (e.g. if protein is fixed)
%     LabelName      name of spin label, 'R1' (default) or 'TOAC'
%     HMM            structure, output from 'mdhmm'
%      .TransProb    transition probability matrix
%      .eqDistr      equilibrium distribution vector
%      .mu           center vectors of states
%      .Sigma        covariance matrices of states
%      .viterbiTraj  Viterbi state trajectory (most likely given the dihedrals)
%      .tauRelax     relaxation times of HMM
%      .logLik       log-likelihood of HMM during optimization
%
%   Output:
%     B              magnetic field vector (mT)
%     spc            EPR spectrum
%     TDSignal       FID time-domain signal
%     t              time axis (in s)


function varargout = cardamom(Sys,Exp,Par,Opt,MD)

if nargin==0, help(mfilename); return; end

% Check expiry date
error(eschecker);

% Check Matlab version
error(chkmlver);

% Preprocessing
% -------------------------------------------------------------------------
if nargin<3
  error('At least 3 inputs (Sys,Exp,Par) are required.');
end
if nargin>5
  error('At most 5 inputs (Sys,Exp,Par,Opt,MD) are possible.');
end
if nargin<4, Opt = struct; end
if nargin<5, MD = []; end

switch nargout
  case 0 % plotting
  case 1 % spc
  case 2 % B,spc
  case 3 % B,spc,TDSignal
  case 4 % B,spc,TDSignal,t
  otherwise
    error('Incorrect number of output arguments.');
end

if ~isfield(Opt,'Verbosity')
  Opt.Verbosity = 0; % Log level
end


global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

% Check Sys
% -------------------------------------------------------------------------

[Sys,err] = validatespinsys(Sys);
error(err);

if Sys.nElectrons>1, error('cardamom does not support more than one electron.'); end
if Sys.nNuclei>1, error('cardamom does not support more than one nucleus.'); end

if isfield(Sys, 'lw')
  isBroadening = any(Sys.lw>0);
else
  isBroadening = false;
end

% Check Exp
%-------------------------------------------------------------------------------

logmsg(1,'Experimental settings:');
[Exp,FieldSweep,CenterField,CenterFreq,Sweep] = validate_exp('cardamom',Sys,Exp);

if FieldSweep
  omega0 = 2*pi*Exp.mwFreq*1e9;  % GHz -> rad s^-1
else
  omega0 = 2*pi*CenterFreq*1e9;  % GHz -> rad s^-1
end

% Check local dynamics models
%-------------------------------------------------------------------------------
if ~isfield(Par,'Model')
  error('Please specify a simulation model in Par.Model.')
end

if ~ischar(Par.Model)
  error('Par.Model must be a string.')
end

switch Par.Model
  case {'diffusion','jump'}
    if ~isempty(MD)
      error('For stochastic diffusion/jump model, do not provide an MD trajectory.');
    end
    
  case {'MD-direct','MD-HBD','MD-HMM'}
    if isempty(MD)
      error('For the model ''%s'', MD simulation information must be provided in the input argument ''MD''.',Par.Model);
    end
  otherwise
    error('Setting ''%s'' for Par.Model not recognized.', Par.Model);
end
LocalDynamicsModel = Par.Model;
useMDdirect = strcmp(LocalDynamicsModel,'MD-direct');
useMD = ~isempty(MD);

% Check MD
%-------------------------------------------------------------------------------
if useMD
  logmsg(1,'  using MD trajectory data');
  
  if ~isfield(MD,'dt')
    error('The MD trajectory time step MD.dt must be given.')
  end
    
  if ~isfield(MD,'DiffGlobal')
    MD.DiffGlobal = [];
  end
  
  if ~isfield(MD,'LabelName')
    MD.LabelName = 'R1';
  end
  
  if ~isfield(MD,'FrameTraj')
    error('The spin label frame trajectory MD.FrameTraj must be given.');
  end
  if size(MD.FrameTraj,1)~=3 || size(MD.FrameTraj,2)~=3
    error('Frame trajectory in MD must be of size (3,3,nSteps,nTraj).');
  end
  
  % Swap last two dimensions if size (...,nTraj,nSteps) is given, to get
  % size (...,nSteps,nTraj)
  dimsSwapped = size(MD.FrameTraj,3)<size(MD.FrameTraj,4);
  if dimsSwapped
    perm34 = @(x)permute(x,[1 2 4 3]);
    MD.FrameTraj = perm34(MD.FrameTraj);
    MD.FrameTrajwrtProt = perm34(MD.FrameTrajwrtProt);
    MD.dihedrals = perm34(MD.dihedrals);
  end
  
  if ~isfield(MD,'FrameTrajwrtProt')
    error('The spin label frame trajectory MD.FrameTrajwrtProt must be given.');
  end
  
  if ~isfield(MD,'removeGlobal')
    MD.removeGlobal = true;
  end
  
  if MD.removeGlobal
    MD.RTraj = MD.FrameTrajwrtProt;
  else
    MD.RTraj = MD.FrameTraj;
  end
  
  MD.nTraj = size(MD.RTraj,4);  % number of trajectories
  MD.nSteps = size(MD.RTraj,3); % number of time steps
  if MD.nTraj~=1
    error('Can only handle MD data with a single trajectory.');
  end
    
  % Check for orthogonality of rotation matrices
  RTrajInv = permute(MD.RTraj,[2,1,3,4]);
  if ~allclose(multimatmult(MD.RTraj,RTrajInv),...
               repmat(eye(3),1,1,size(MD.RTraj,3),size(MD.RTraj,4)),...
               1e-10)
    error('Rotation matrices in frame trajectory are not orthogonal.')
  end
  clear RTrajInv
  
  logmsg(1,'    label: %s',MD.LabelName);
  logmsg(1,'    number of trajectories: %d',MD.nTraj);
  logmsg(1,'    number of time steps: %d',MD.nSteps);
  logmsg(1,'    size of time step: %g ps',MD.dt/1e-12);
  logmsg(1,'    remove global diffusion: %d',MD.removeGlobal);
  
  % Build Markov state model
  if strcmp(LocalDynamicsModel,'MD-HMM')
    logmsg(1,'Building HMM model');
    if isfield(MD,'HMM')
      logmsg(1,'  using provided HMM parameters');
      
      HMM = MD.HMM;
      if ~isfield(HMM,'TransProb')
        error('The transition probability matrix must be provided in HMM.transmat.')
      end
      if ~isfield(HMM,'nLag')
        error('The time lag (an integer multiple of the MD time step) must be provided in HMM.nLag.')
      end
      if ~isfield(HMM,'eqDistr')
        error('The equilibrium state distribution must be provided in HMM.eqdistr.')
      end
      if ~isfield(HMM,'viterbiTraj')
        error('The Viterbi trajectory must be provided in HMM.viterbiTraj.')
      end
      
    else
      logmsg(1,'  constructing HMM from MD');

      nLag = round(MD.tLag/MD.dt);
      HMM = mdhmm(MD.dihedrals,MD.dt,MD.nStates,nLag,Opt);
      
    end
    % provides HMM.transmat, HMM.eqdistr, HMM.viterbiTraj, etc
    MD.viterbiTraj = HMM.viterbiTraj.';
    MD.nStates = HMM.nStates;
    
    % Set the Markov chain time step based on the (scaled) sampling lag time
    Par.dtSpatial = MD.tLag;
  end
  
  % Estimate rotational diffusion time constant (used in the density propagation)
  FrameAcorr = squeeze(autocorrfft(squeeze(MD.FrameTraj.^2), 3, 1, 1, 1));
  N = round(MD.nSteps/2);
  time = linspace(0, N*MD.dt, N);
  tauR = max(cumtrapz(time,FrameAcorr(:,1:N),2),[],2);
  tauR = mean(tauR);
  DiffLocal = 1/6/tauR;
  MD.tauR = tauR;  
  
end


% Check dynamics and ordering
%-------------------------------------------------------------------------------

isDiffSim = (useMD && strcmp(LocalDynamicsModel,'MD-HBD')) || ...
   strcmp(LocalDynamicsModel,'diffusion');

dynamInfoGiven = ( isfield(Sys,'tcorr') || isfield(Sys,'Diff') ...
                   || isfield(Sys,'logtcorr') || isfield(Sys,'logDiff') );
if useMD && strcmp(LocalDynamicsModel,'MD-HBD') && ~dynamInfoGiven
  % estimate rotational diffusion tensor
  % currently only supports a single MD trajectory
  % TODO: make this work for multiple MD trajectories
  logmsg(1,'Sys.Diff not specified, estimating from MD trajectory data')
  stopFitT = floor(MD.nSteps/2)*MD.dt;
  [Sys.Diff, ~, ~] = runprivate('cardamom_estimatedifftensor',...
                                squeeze(MD.FrameTraj), MD.dt, stopFitT);
  logmsg(1,'Estimated Sys.Diff eigenvalues:  (%g, %g, %g) rad^2/us',Sys.Diff/1e6);
end

Dynamics = validate_dynord('cardamom',Sys,FieldSweep,isDiffSim);

if isDiffSim
  DiffLocal = Dynamics.Diff;
  Sys.Diff = DiffLocal;
end

if useMD
  Dynamics.DiffGlobal = MD.DiffGlobal;
end
includeGlobalDynamics = ~isempty(Dynamics.DiffGlobal);

if isfield(Sys,'Potential')
  useLocalPotential = true;
  LocalPotential = Sys.Potential;
else
  useLocalPotential = false;
end

% Check Opt
%-------------------------------------------------------------------------------

fastMethodApplicable = Sys.nElectrons==1 && Sys.S==1/2 && Sys.nNuclei<=1;
if ~isfield(Opt,'Method')
  if fastMethodApplicable
    Opt.Method = 'fast';
  else
    Opt.Method = 'ISTOs';
  end
end
if strcmp(Opt.Method,'fast')
  if ~fastMethodApplicable
    error('Opt.Method = ''fast'' is not applicable for this spin system.');
  end
end

if ~isfield(Opt,'FFTWindow')
  Opt.FFTWindow = true;
end
FFTWindow = Opt.FFTWindow;

if ~isfield(Opt,'specCon')
  Opt.specCon = false;
end
checkConvergence = Opt.specCon;

% Lag time for sliding window processing (used in MD-direct only)
if ~isfield(Opt,'LagTime')
  Opt.LagTime = 2e-9; % seconds
end

% Check Par
%-------------------------------------------------------------------------------

% Require Par.dtSpin
if ~isfield(Par,'dtSpin')
  error('Par.dtSpin (spin propagation time step) must be given.');
end

% Require and check Par.nSteps
if ~isfield(Par,'nSteps')
  error('Par.nSteps must be given.');
end
if ~isnumeric(Par.nSteps) || numel(Par.nSteps)~=1 || ~isreal(Par.nSteps) || ...
    mod(Par.nSteps,1) || Par.nSteps<1
  error('Par.nSteps must be a positive integer.');
end

% Check Par.dtSpatial
if isfield(Par,'dtSpatial')
  if useMDdirect
    error('For MD-direct simulations, Par.dtSpatial is not allowed. The time step is taken from the MD input.');
  end
  if Par.dtSpin<Par.dtSpatial
    error('The spatial dynamics time step Par.dtSpatial must be less than or equal to the spin dynamics time step Par.dtSpin.')
  end
else
  if ~useMDdirect
    error('The spatial dynamics time step Par.dtSpatial must be specified when using an %s model.',LocalDynamicsModel);
  end
  Par.dtSpatial = MD.dt;
end

% Make sure Par.dtSpin is an (approx.) integer multiple of Par.dtSpatial
r = Par.dtSpin/Par.dtSpatial;
if abs(r-round(r))>1e-4 || r<1
  error('The spin propagation time step (Par.dtSpin) must be a multiple of the trajectory time step (Par.dtSpatial or MD.dt).');
end
r = round(r);
Par.dtSpatial = Par.dtSpin/r;
Par.BlockLength = r; % used for block averaging

nStepsSpin = Par.nSteps;
nStepsSpatial = nStepsSpin*round(Par.dtSpin/Par.dtSpatial);
dtSpin = Par.dtSpin;
dtSpatial = Par.dtSpatial;

% Set default number of (stochastic) trajectories
if ~useMDdirect && ~isfield(Par,'nTraj')
  Par.nTraj = 100; 
end

% Determine whether to process single long MD trajectory into multiple short
% MD trajectories
if useMDdirect
  Par.lag = ceil(Opt.LagTime/Par.dtSpin);
  nBlocks = floor(MD.nSteps/Par.BlockLength);
  if nBlocks < Par.nSteps
    error('MD trajectory is too short (%g ns) for the required FID length (%g ns.',...
      MD.nSteps*MD.dt,Par.nSteps*Par.dtSpin);
  end
  Par.nTraj = floor((nBlocks-Par.nSteps)/Par.lag) + 1;
else
  Par.lag = 0;
end

% Set default number of orientations
if ~isfield(Par,'nOrients')
  Par.nOrients = 100;
end

logmsg(1,'Parameter settings:');
logmsg(1,'  Local dynamics model:   ''%s''',LocalDynamicsModel);
if includeGlobalDynamics
  logmsg(1,'  Global correlation time:  %g rad^2/us',Dynamics.DiffGlobal/1e6);
else
  logmsg(1,'  Global correlation time:  none');
end  
logmsg(1,'  Number of orientations: %d',Par.nOrients);
logmsg(1,'  Number of trajectories: %d',Par.nTraj);
logmsg(1,'  Spin propagation:       %d steps of %g ns (%g ns total)',nStepsSpin,dtSpin/1e-9,nStepsSpin*dtSpin/1e-9);
logmsg(1,'  Spatial propagation:    %d steps of %g ns (%g ns total)',nStepsSpatial,dtSpatial/1e-9,nStepsSpatial*dtSpatial/1e-9);
logmsg(1,'  Lag time:               %g MD steps',Par.lag);


% Check local dynamics model
%-------------------------------------------------------------------------------
switch LocalDynamicsModel
  case 'diffusion'
    
    DiffLocal = Dynamics.Diff;
    
  case 'jump'
    
  case {'MD-direct','MD-HBD','MD-HMM'} % TODO process RTraj based on size of input
    
    if ~isfield(Par,'nOrients')
      error('Par.nOrients must be specified for an MD-based model.')
    end
    
    RTrajLocal = MD.RTraj;
    qTrajLocal = rotmat2quat(RTrajLocal);
    
    switch LocalDynamicsModel
      case 'MD-HBD'
        
        % calculate orienting potential energy function
        qTemp = squeeze(rotmat2quat(MD.FrameTraj));
        [phi,theta,psi] = quat2euler(qTemp,'active');
        clear qTemp
        
        phi = phi + 2*pi*(phi<0);
        psi = psi + 2*pi*(psi<0);
        
        nBins = 90;
        phiEdges = linspace(0, 2*pi, nBins+1);
        thetaEdges = linspace(0, pi, nBins/2+1);
        psiEdges = linspace(0, 2*pi, nBins+1);
        
        pdf = histcountsn([phi(:),theta(:),psi(:)],{phiEdges,thetaEdges,psiEdges});
        pdf = smooth3(pdf,'gaussian');
        pdf(pdf<1e-14) = 1e-14;  % put a finite floor on histogram
        useLocalPotential = true;
        LocalPotential = -log(pdf);
        
      case 'MD-HMM'
        
        offset = 1;
        RTrajLocal = RTrajLocal(:,:,offset:HMM.nLag:end,:);
        qTrajLocal = rotmat2quat(RTrajLocal);
        
    end
    
    % Delete very large arrays that are no longer needed
    fields = {'FrameTraj','RTraj','FrameTrajwrtProt','RProtDiff'};
    for iField = 1:numel(fields)
      if isfield(MD,fields{iField})
        MD = rmfield(MD,fields{iField});
      end
    end
    
end

% Generate grids for powder averaging in the lab frame
%-------------------------------------------------------------------------------

nOrientations = Par.nOrients;
Orientations = [];
if isfield(Par,'Orients')
  Orientations = Par.Orients;
  if isvector(Orientations)
    Orientations = repmat(Orientations(:),[1,nOrientations]);
  end
end

% Set up orientational grid for powder averaging
if ~isempty(Orientations)
  gridPhi = Orientations(1,:);
  gridTheta = Orientations(2,:);
  weight = ones(size(gridPhi));
  weight = weight/sum(weight);
else
  if checkConvergence
    % Generate Sobol sequence over a spiral grid
    nOrientations = ceil(nOrientations/2);
    skip = 0;
    gridPts = 2*cardamom_sobol_generate(1,nOrientations,skip)-1;
    gridPhi = sqrt(pi*nOrientations)*asin(gridPts);
    gridTheta = acos(gridPts);
    weight = ones(size(gridPhi));
  else
    % Generate a spiral grid over the full sphere
    gridPts = linspace(-1,1,nOrientations);
    gridPhi = sqrt(pi*nOrientations)*asin(gridPts);
    gridTheta = acos(gridPts);
    weight = ones(size(gridPhi));
    weight = weight/sum(weight);
    % Generate a triangular grid over the upper hemisphere (since EPR spectra
    % are invariant under inversion)
    %GridSize = ceil(sqrt(nOrientations));
    %grid = sphgrid('Ci',GridSize);
    %gridPhi = grid.phi;
    %gridTheta = grid.theta;
    %weight = grid.weights;
    %weight = weight/sum(weight);
  end
end


% Run simulation
%-------------------------------------------------------------------------------
logmsg(1,'Quantum propagation method: ''%s''',Opt.Method);
logmsg(1,'Running simulation');

clear cardamom_propagatedm % to clear persistent variables in function

iter = 1;
spcArray = [];
nOrientsTot = 0;

t = linspace(0, nStepsSpin*Par.dtSpin, nStepsSpin).';

converged = false;
spcLast = 0;
while ~converged
  tic
  
  % store signals from trajectories in cell array
  TDSignal = {};
  tCell = [];
  
  % temporary cells to store intermediate results
  iTDSignal = cell(1,nOrientations);
  itCell = cell(1,nOrientations);
  
  updateuser(0);
  for iOri = 1:nOrientations
    
    % Generate trajectory of local dynamics
    switch LocalDynamicsModel
      
      case 'diffusion'
        
        Sys.Diff = Dynamics.Diff;
        Par.dt = dtSpatial;
        Par.nSteps = nStepsSpatial;
        if useLocalPotential
          Sys.Potential = LocalPotential;
          Par.nSteps = 2*nStepsSpatial;
        end
        [~, RTrajLocal, qTrajLocal] = stochtraj_diffusion(Sys,Par,Opt);
        if useLocalPotential
          RTrajLocal = RTrajLocal(:,:,nStepsSpatial+1:end,:);
          qTrajLocal = qTrajLocal(:,:,nStepsSpatial+1:end,:);
        end
        
      case 'jump'
        
        Par.dt = dtSpatial;
        Par.nSteps = nStepsSpatial;
        [~, RTrajLocal, qTrajLocal] = stochtraj_jump(Sys,Par,Opt);
        
      case 'MD-direct'
        
        % the MD trajectories are not changing with orientation,
        % RTraj and qTraj were processed earlier outside of the
        % orientation loop
        
      case 'MD-HBD'
        
        Sys.Diff = DiffLocal;
        Par.nSteps = nStepsSpatial;
        if useLocalPotential
          Sys.Potential = LocalPotential;
          Par.nSteps = 2*nStepsSpatial;
        end
        Par.dt = dtSpatial;
        [~, RTrajLocal, qTrajLocal] = stochtraj_diffusion(Sys,Par,Opt);
        if useLocalPotential
          RTrajLocal = RTrajLocal(:,:,nStepsSpatial+1:end,:);
          qTrajLocal = qTrajLocal(:,nStepsSpatial+1:end,:);
        end
        
      case 'MD-HMM'
        
        Sys.TransProb = HMM.TransProb;
        Par.dt = dtSpatial;
        Par.nSteps = 2*nStepsSpatial;
        CumulDist = cumsum(HMM.eqDistr)/sum(HMM.eqDistr);
        for k = 1:Par.nTraj
          Par.StatesStart(k) = find(CumulDist>rand(),1);
        end
        Opt.statesOnly = true;
        [~, stateTraj] = stochtraj_jump(Sys,Par,Opt);
        stateTraj = stateTraj(nStepsSpatial+1:end,:);

        Par.stateTraj = stateTraj;
        
    end
    
    if strcmp(Opt.Method,'ISTOs')
      Par.qTraj = qTrajLocal;
    else
      Par.RTraj = RTrajLocal;
    end
    
    qLab = euler2quat(0, gridTheta(iOri), gridPhi(iOri), 'active');
    qLab = repmat(qLab,[1,nStepsSpin,Par.nTraj]);
    
    % Generate trajectory of global dynamics with a time step equal to that 
    % of the spin propagation (these rotations will be performed AFTER
    % the time-dependent interaction tensors are calculated and possibly 
    % averaged)
    if includeGlobalDynamics
      Sys_.Diff = Dynamics.DiffGlobal;
      Par.dt = dtSpin;
      Par.nSteps = nStepsSpin;
      [~, ~, qTrajGlobal] = stochtraj_diffusion(Sys_,Par,Opt);
      % Combine global trajectories with starting orientations
      qLab = quatmult(qLab,qTrajGlobal);
    end
    
    if strcmp(Opt.Method,'ISTOs')
      Par.qLab = qLab;
    else
      Par.RLab = quat2rotmat(qLab);
    end
    
    % Propagate the density matrix
    Par.nSteps = nStepsSpin;
    Par.dtSpin = dtSpin;
    Par.dtSpatial = dtSpatial;
    Sprho = cardamom_propagatedm(Sys,Par,Opt,MD,omega0,CenterField);
    
    % Calculate the time-domain signal, i.e. the expectation value of S_{+}
    iTDSignal{1,iOri} = 0;
    for k = 1:size(Sprho,1)
      iTDSignal{1,iOri} = iTDSignal{1,iOri} + squeeze(Sprho(k,k,:));
    end
    iTDSignal{1,iOri} = weight(iOri)*iTDSignal{1,iOri};
    
    if Opt.Verbosity
      updateuser(iOri,nOrientations,true);
    end
    
    itCell{1,iOri} = t;
    
  end
  
  % Store simulations at new starting orientations from each iteration
  TDSignal = cat(2, TDSignal, iTDSignal);
  tCell = cat(2, tCell, itCell);

  % Perform FFT
  % -------------------------------------------------------------------------

  % windowing
  if FFTWindow
    hamming = cellfun(@(x) 0.54 + 0.46*cos(pi*x/max(x)), tCell, 'UniformOutput', false);
    hann = cellfun(@(x) 0.5*(1 + cos(pi*x/max(x))), tCell, 'UniformOutput', false);
    Window = hann;
    TDSignal = cellfun(@times, TDSignal, Window, 'UniformOutput', false);
  end

  % zero padding for FFT to ensure sufficient B-field resolution (at most 0.1 G)
  Bres = 0.1; % G
  tReq = 1/(mt2mhz(Bres/10)*1e6); % mT -> s

  tMax = max(cellfun(@(x) max(x), tCell));
  if tMax<tReq
    M = ceil(tReq/Par.dtSpin);
  else
    M = ceil(tMax/Par.dtSpin);
  end
  TDSignal = cell2mat(cellfun(@(x) zeropad(x, M), TDSignal, 'UniformOutput', false));
  tLong = linspace(0, M*Par.dtSpin, M).';

  % convolute for linewidth broadening
  if isBroadening
    if Sys.lw(1)>0
      % Gaussian broadening
      fwhm = mt2mhz(Sys.lw(1))*1e6;  % mT -> Hz
      alpha = pi^2*fwhm^2/(4*log(2));
      TDSignal = bsxfun(@times,exp(-alpha*tLong.^2),TDSignal);
    end
    if numel(Sys.lw)==2 && Sys.lw(2)>0
      % Lorentzian broadening
      TL = Dynamics.T2;
      TDSignal = bsxfun(@times,exp(-tLong/TL),TDSignal);
    end
  end
  
  % Multiply by t for differentiation and take the FFT
  spcArray = cat(2, spcArray, imag(fftshift(fft(bsxfun(@times,TDSignal,tLong), [], 1))));
  spcNew = mean(spcArray,2);
  
  if checkConvergence
    span = max(spcNew)-min(spcNew);
    rmsdNew = sqrt(mean((spcNew-spcLast).^2))/span;
    if rmsdNew<2e-2
      converged = true;
    elseif iter>2
      rmsdRelativeChange = abs((rmsdNew-rmsdLast)/rmsdLast);
      converged = rmsdRelativeChange<0.5;
    end
  else
    converged = true;
  end
  
  if ~converged
    % double the number of orientations in powder averaging
    msg = sprintf('Convergence not achieved. Propagation is being extended.\n');
    if Opt.Verbosity
      fprintf(msg);
    end
    rmsdLast = rmsdNew;
    spcLast = spcNew;  % store for comparison after next iteration completes
    nOrientsTot = nOrientsTot + nOrientations;
    skip = iter*nOrientsTot;  % seed Sobol sequence generator for next iteration
    
    nOrientations = nOrientsTot;
    if ~isempty(Orientations)
      gridPhi = repmat(Orientations(:,1),[2^(iter-1),1]);
      gridTheta = repmat(Orientations(:,2),[2^(iter-1),1]);
    else
      gridPts = 2*cardamom_sobol_generate(1,nOrientations,skip)-1;
      gridPhi = sqrt(pi*nOrientations)*asin(gridPts);
      gridTheta = acos(gridPts);
    end
    iter = iter + 1;
  else
    spcAvg = spcNew;
    minsTot = floor(toc/60);
    if Opt.Verbosity
      msg = sprintf('Done!\nTotal simulation time: %0d:%02.0f\n',minsTot,mod(toc,60));
     fprintf(msg);
    end
  end

end

clear RTraj qTraj
% Par = rmfield(Par,{'RTraj','qTraj'});

if strcmp(LocalDynamicsModel(1:2), 'MD')
  % these variables can take up a lot of memory and might prevent the user 
  % from implementing a fine enough grid for powder averaging 
  clear MD
end

freq = 1/(Par.dtSpin*M)*(-M/2:M/2-1);  % TODO check for consistency between FieldSweep and FFT window/resolution

% center the spectrum around the isotropic component of the g-tensor
if FieldSweep
  fftAxis = mhz2mt(freq/1e6+Exp.mwFreq*1e3, mean(Sys.g));  % MHz -> mT, note use of g0, not ge
else
  fftAxis = freq/1e9+mt2mhz(Exp.Field,mean(Sys.g))/1e3;  % GHz
end

% set up horizontal sweep axis
if FieldSweep
  xAxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);  % field axis, mT
else
  xAxis = linspace(Exp.mwRange(1),Exp.mwRange(2),Exp.nPoints);  % field axis, GHz
end

% interpolate over horizontal sweep range
if FieldSweep
  outspc = interp1(fftAxis, spcAvg, xAxis);
else
  spcAvg = spcAvg(end:-1:1);   % reverse the axis for frequency sweep
  outspc = interp1(fftAxis,spcAvg,xAxis,'spline',0);
  outspc = cumtrapz(xAxis(end:-1:1),outspc);  % frequency sweeps outputs the absorption
end

% average over trajectories for time-domain signal output
TDSignal = mean(TDSignal, 2);
TDSignal = TDSignal(1:Par.nSteps);

% Plotting, output
%-------------------------------------------------------------------------------
switch nargout
case 0
  cla
  if FieldSweep
    if xAxis(end)<10000
      plot(xAxis,outspc);
      xlabel('magnetic field (mT)');
    else
      plot(xAxis/1e3,outspc);
      xlabel('magnetic field (T)');
    end
    axis tight
    ylabel('intensity (arb.u.)');
    title(sprintf('%0.8g GHz',Exp.mwFreq));
  else
    if xAxis(end)<1
      plot(xAxis*1e3,outspc);
      xlabel('frequency (MHz)');
    else
      plot(xAxis,outspc);
      xlabel('frequency (GHz)');
    end
    axis tight
    ylabel('intensity (arb.u.)');
    title(sprintf('%0.8g mT',Exp.Field));
  end
case 1
  varargout = {outspc};
case 2
  varargout = {xAxis,outspc};
case 3
  varargout = {xAxis,outspc,TDSignal};
case 4
  varargout = {xAxis,outspc,TDSignal,t};
end

clear global EasySpinLogLevel
    
end
%===============================================================================


% Helper functions
% -------------------------------------------------------------------------

function updateuser(iOrient,nOrient,reverse)
% Update user on progress

persistent reverseStr

if iOrient==0
  reverseStr = '';
  return
end

secsElapsed = toc;
minsElapsed =  floor(secsElapsed/60);
avgTime = secsElapsed/iOrient;
secsLeft = (nOrient - iOrient)*avgTime;
minsLeft = floor(secsLeft/60);

msg2 = sprintf('  Orientation %d/%d:\n',iOrient,nOrient); 
msg3 = sprintf('   Time elapsed %02d:%02d:%02.0f (%0.3g s/orientation)\n', floor(minsElapsed/60), mod(minsElapsed,60), mod(secsElapsed,60),avgTime);
msg4 = sprintf('   remaining %02d:%02d:%02.0f\n', floor(minsLeft/60), mod(minsLeft,60), mod(secsLeft,60));
msg = [msg2 msg3 msg4];

if reverse
  fprintf([reverseStr msg]);
else
  fprintf(msg);
end

reverseStr = repmat(sprintf('\b'), 1, length(msg));

end

function  y = zeropad(x, M)
  N = length(x);
  if iscolumn(x), y = [x; zeros(M-N, 1)]; end
  if isrow(x), y = [x, zeros(1, M-N)]; end
end
