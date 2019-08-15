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
%                       of the orientational potential for that orientation. The
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
%     Model          string 
%                    the model for spin label dynamics
%                    'diffusion'
%                    'jump'
%                    'MD-direct': use molecular orientations in MD
%                      trajectories directly as input for simulating the
%                      spectrum
%                    'MD-HBD': coarse grain the MD trajectories by using 
%                      the Euler angle probability distribution (for 
%                      pseudopotential) from the spin label's orientations 
%                      to perform further stochastic rotational dynamics 
%                      simulations
%                    'MD-HMM': coarse grain the MD trajectories by using 
%                      the spin label's side chain dihedral angles to build 
%                      a hidden Markov model model to perform further 
%                      stochastic jump dynamics simulations
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
%                    fast: propagate the density matrix using an 
%                      analytical expression for the matrix exponential in 
%                      the m_S=-1/2 subspace (S=1/2 with up to one nucleus)
%                    ISTOs: propagate the density matrix using
%                      irreducible spherical tensor operators (general, slower)
%
%     FFTWindow       1: use a Hamming window (default), 0: no window
%
%     nTrials        integer
%                    number of initialization trials for k-means
%                    clustering; used for the Markov method
% 
%     LagTime        lag time for sliding window processing
%
%
%   MD: structure with molecular dynamics simulation parameters
%
%     dt             double
%                    time step (in s) for saving MD trajectory snapshots
%
%     tLag           double
%                    time lag (in s) for sampling the MD trajectory to 
%                    determine states and transitions, used for the hidden
%                    Markov model
%
%     nStates        number of states in the hidden Markov model
%
%     DiffGlobal     double (optional)
%                    Diffusion coefficient for isotropic global rotational
%                    diffusion (s^-1)
% 
%     removeGlobal   integer
%                    1: (default) remove protein global diffusion
%                    0: no removal (e.g. if protein is fixed)
% 
%     LabelName      name of spin label, 'R1' (default) or 'TOAC'
%
%     HMM            structure, output from 'mdhmm'
%      .TransProb    transition probability matrix
%      .eqDistr      equilibrium distribution vector
%      .mu           center vectors of states
%      .Sigma        covariance matrices of states
%      .viterbiTraj  Viterbi state trajectory (most likely given the dihedrals)
%      .tauRelax     relaxation times of HMM
%      .logLik       log-likelihood of HMM during optimization
%                    
%
%   Output:
%     B              numeric, size = (2*nSteps,1) 
%                    magnetic field (mT)
%
%     spc            numeric, size = (2*nSteps,1)
%                    derivative EPR spectrum
%
%     TDSignal       numeric, size = (2*nSteps,1)
%                    time-domain signal, <S_+>
%
%     t              numeric, size = (2*nSteps,1)
%                    simulation time axis (in s)

%    Opt.ExpMethod   method for computing matrix exponential

function varargout = cardamom(Sys,Exp,Par,Opt,MD)

% Preprocessing
% -------------------------------------------------------------------------

switch nargin
  case 0
    help(mfilename); return;
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
  case 3 % B,spc,TDSignal
  case 4 % B,spc,TDSignal,t
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

% Check local dynamics models
%-------------------------------------------------------------------------------
if ~isfield(Par,'Model')
  error('Please specify a simulation model using Par.Model.')
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
    error('The spin label frame trajectory MD.FrameTraj must be given.')
  end
  
  if ~isfield(MD,'removeGlobal')
    MD.removeGlobal = true;
  end  
  if MD.removeGlobal
    MD.RTraj = MD.FrameTrajwrtProt;
  else
    MD.RTraj = MD.FrameTraj;
  end
  
  if size(MD.RTraj,1)~=3 || size(MD.RTraj,2)~=3
    error('Frame trajectory in MD must be of size (3,3,nTraj,nSteps).');
  end
  MD.nTraj = size(MD.RTraj,3);  % number of trajectories; assumed to be one for now
  MD.nSteps = size(MD.RTraj,4); % number of time steps
    
  % Check for orthogonality of rotation matrices
  RTrajInv = permute(MD.RTraj,[2,1,3,4]);
  if ~allclose(multimatmult(MD.RTraj,RTrajInv),...
               repmat(eye(3),1,1,size(MD.RTraj,3),size(MD.RTraj,4)),...
               1e-10)
    error('Rotation matrices in frame trajectory are not orthogonal.')
  end
  clear RTrajInv
  
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

      nLag = MD.tLag/MD.dt;
      HMM = mdhmm(MD.dihedrals,MD.dt,MD.nStates,nLag,Opt);
      
    end
    % provides HMM.transmat, HMM.eqdistr, HMM.viterbiTraj, etc
    MD.viterbiTraj = HMM.viterbiTraj.';
    MD.nStates = HMM.nStates;

    % Set the Markov chain time step based on the (scaled) sampling lag time
    Par.dt = MD.tLag;
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

if useMD
  isDiffSim = strcmp(LocalDynamicsModel,'MD-HBD');
elseif strcmp(LocalDynamicsModel,'jump')
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
  if isfield(MD, 'Potential')
    isGlobalPotential = true;
    GlobalPotential = MD.Potential;  % TODO fully implement this in Dynamics-checking and documentation
  else
    isGlobalPotential = false;
  end
end

if isfield(Sys,'Potential')
  isLocalPotential = true;
  LocalPotential = Sys.Potential;
else
  isLocalPotential = false;
end

% Check Opt
%-------------------------------------------------------------------------------

fastMethodApplicable = Sys.nElectrons==1 && Sys.nNuclei<=1;
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

if ~isfield(Opt,'ExpMethod')
  Opt.ExpMethod = 'eig';
end

% Lag time for sliding window processing (used in MD-direct only)
if ~isfield(Opt,'LagTime')
  Opt.LagTime = 2e-9; % seconds
end

% Check Par
%-------------------------------------------------------------------------------

% Set default number of (stochastic) trajectories
if ~isfield(Par,'nTraj') && ~strcmp(LocalDynamicsModel,'MD-direct')
  Par.nTraj = 100; 
end

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
  elseif strcmp(LocalDynamicsModel,'jump') && ~isfield(Par,'dt')
    error('The time step Par.dt must be specified when using an jump model.')
  elseif strcmp(LocalDynamicsModel,'MD') && ~isfield(Par,'dt')
    error('The time step Par.dt must be specified when using an MD model.')
  end
  Par.Dt = Par.dt;
  if isfield(Par,'nSteps')
    nSteps = Par.nSteps;
  else
    nSteps = ceil(250e-9/Par.dt);
  end
  nStepsStoch = nSteps;
  nStepsQuant = nSteps;
end

dtQuant = Par.Dt;
dtStoch = Par.dt;

if useMDdirect
  if Par.Dt<MD.dt
    error('Quantum time step (Par.Dt) cannot be smaller than MD time step (MD.dt).');
  end
end

% Calculate block length for block averaging (done if block length > 1)
if useMDdirect
  traj_dt = MD.dt;
else
  traj_dt = Par.dt;
end
Par.BlockLength = ceil(Par.Dt/traj_dt);

% Determine whether to process single long MD trajectory into multiple short
% MD trajectories
if useMDdirect
  Par.lag = ceil(Opt.LagTime/Par.Dt);
  nBlocks = floor(MD.nSteps/Par.BlockLength);
  if Par.nSteps<nBlocks
    % Par.nSteps not changed from user input
    Par.nTraj = floor((nBlocks-Par.nSteps)/Par.lag) + 1;
  else
    Par.nSteps = nBlocks;
    Par.nTraj = 1;
  end
end

% Set default number of orientations
if ~isfield(Par,'nOrients')
  Par.nOrients = 100;
end

logmsg(1,'Parameter settings:');
logmsg(1,'  Local dynamics model:   ''%s''',LocalDynamicsModel);
logmsg(1,'  Number of trajectories: %d',Par.nTraj);
logmsg(1,'  Number of orientations: %d',Par.nOrients);
logmsg(1,'  Quantum propagation:    %d steps of %g ns',nStepsQuant,dtQuant/1e-9);
logmsg(1,'  Spatial propagation:    %d steps of %g ns',nStepsStoch,dtStoch/1e-9);


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
        phiBins = linspace(0, 2*pi, nBins+1);
        thetaBins = linspace(0, pi, nBins/2+1);
        psiBins = linspace(0, 2*pi, nBins+1);
        
        pdf = histcnd([phi(:),theta(:),psi(:)],{phiBins,thetaBins,psiBins});
        
        pdf(end,:,:) = pdf(1,:,:);  % FIXME why does it truncate to zero in the phi direction?
        pdf = smooth3(pdf,'gaussian');
        pdf(pdf<1e-14) = 1e-14;  % put a finite floor on histogram
        isLocalPotential = true;
        LocalPotential = -log(pdf);
        
      case 'MD-HMM'
        
        RTrajLocal = RTrajLocal(:,:,:,1:HMM.nLag:end);
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
Orients = [];
if isfield(Par,'Orients')
  Orients = Par.Orients;
  if isvector(Orients)
    % assure that Orients is a column vector
    Orients = Orients(:);
    if nOrientations>1
      % if only one orientation, then repeat it nOrients times
      Orients = repmat(Orients,[1,nOrientations]);
    end
  end
end

% Set up orientational grid for powder averaging
if ~isempty(Orients)
  gridPhi = Orients(1,:);
  gridTheta = Orients(2,:);
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
    %nKnots = ceil(sqrt(nOrientations));
    %[gridPhi,gridTheta,weight] = sphgrid('Ci',nKnots);
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

t = linspace(0, nStepsQuant*Par.Dt, nStepsQuant);

converged = false;
while ~converged
  tic

  % trajectories might differ in length, so we need cells for allocation
  TDSignal = {};
  tCell = [];

  % temporary cells to store intermediate results
  iTDSignal = cell(1,nOrientations);
  itCell = cell(1,nOrientations);

  updateuser(0);
  for iOrient = 1:nOrientations
    logmsg(1,' Orientation %d/%d',iOrient,nOrientations);
    % Generate trajectory of local dynamics
    switch LocalDynamicsModel
      
      case 'diffusion'
        
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
        
      case 'MD-direct'
            
        % the MD trajectories are not changing, so RTraj and qTraj were
        % processed earlier outside of the orientation loop
            
      case 'MD-HBD'
            
        Sys.Diff = DiffLocal;
        Par.nSteps = nStepsStoch;
%             Sys.Potential = LocalPotential;
        if isLocalPotential
          Sys.Potential = LocalPotential;
          Par.nSteps = 2*nStepsStoch;
        end
        Par.dt = dtStoch;
        [~, RTrajLocal, qTrajLocal] = stochtraj_diffusion(Sys,Par,Opt);
        if isLocalPotential
          RTrajLocal = RTrajLocal(:,:,:,nStepsStoch+1:end);
          qTrajLocal = qTrajLocal(:,:,nStepsStoch+1:end);
        end

      case 'MD-HMM'

        Sys.TransProb = HMM.TransProb;
        Par.dt = dtStoch;
        Par.nSteps = 2*nStepsStoch;
        CumulDist = cumsum(HMM.eqDistr)/sum(HMM.eqDistr);
        for k = 1:Par.nTraj
          Par.StatesStart(k) = find(CumulDist>rand(),1);
        end
        Opt.statesOnly = true;
        [~, stateTraj] = stochtraj_jump(Sys,Par,Opt);
        stateTraj = stateTraj(:,nStepsStoch+1:end);

        Par.stateTraj = stateTraj;
        
    end
    
    if strcmp(Opt.Method,'ISTOs')
      Par.qTraj = qTrajLocal;
    else
      Par.RTraj = RTrajLocal;
    end
    
    % Generate trajectory of global dynamics with a time step equal to that 
    % of the quantum propagation (these rotations will be performed AFTER
    % the time-dependent interaction tensors are calculated and possibly 
    % averaged)
    includeGlobalDynamics = ~isempty(Dynamics.DiffGlobal);
    if includeGlobalDynamics
     Sys.Diff = Dynamics.DiffGlobal;
     if isLocalPotential
       Sys.Potential = [];
     elseif isGlobalPotential
       Sys.Potential = GlobalPotential;
     end
     Par.dt = dtQuant;
     Par.nSteps = nStepsQuant;
     [~, ~, qTrajGlobal] = stochtraj_diffusion(Sys,Par,Opt);
    end
    
    % Combine global trajectories with starting orientations
    qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient), 'active'),...
                  [1,Par.nTraj,nStepsQuant]);
    if includeGlobalDynamics
      qLab = quatmult(qLab,qTrajGlobal);
    end
    
    if strcmp(Opt.Method,'ISTOs')
      Par.qLab = qLab;
    else
      Par.RLab = quat2rotmat(qLab);
    end
    
    % Propagate the density matrix
    Par.nSteps = nStepsQuant;
    Par.Dt = dtQuant;
    Par.dt = dtStoch;
    Sprho = cardamom_propagatedm(Sys,Par,Opt,MD,omega0,CenterField);
    
    % Calculate the time-domain signal, i.e. the expectation value of S_{+}
    iTDSignal{1,iOrient} = 0;
    for k = 1:size(Sprho,1)
      iTDSignal{1,iOrient} = iTDSignal{1,iOrient} + squeeze(Sprho(k,k,:));
    end
    
    iTDSignal{1,iOrient} = weight(iOrient)*iTDSignal{1,iOrient};
    
    if Opt.Verbosity
      updateuser(iOrient,nOrientations,false);
    end
    
    itCell{1,iOrient} = t.';
    
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
    M = ceil(tReq/Par.Dt);
  else
    M = ceil(tMax/Par.Dt);
  end
  TDSignal = cell2mat(cellfun(@(x) zeropad(x, M), TDSignal, 'UniformOutput', false));
  tLong = linspace(0, M*Par.Dt, M).';

  % convolute for linewidth broadening
  if isBroadening
    if Sys.lw(1)>0
      % Gaussian broadening
      fwhm = mt2mhz(Sys.lw(1))*1e6;  % mT -> Hz
      alpha = pi^2*fwhm^2/(4*log(2));
      TDSignal = bsxfun(@times,exp(-alpha*tLong.^2),TDSignal);
    end
    if numel(Sys.lw)==2
      % Lorentzian broadening
      TL = Dynamics.T2; 
      TDSignal = bsxfun(@times,exp(-tLong/TL),TDSignal);
    end
  end
  
  % Multiply by t for differentiation and take the FFT
  spcArray = cat(2, spcArray, imag(fftshift(fft(bsxfun(@times,TDSignal,tLong), [], 1))));
  spcNew = mean(spcArray,2);

  if checkConvergence
    if iter==1
      spcLast = spcNew;
      rmsdNew = NaN;
    else
      span = max(spcNew)-min(spcNew);
      rmsdNew = sqrt(mean((spcNew-spcLast).^2))/span;
      if rmsdNew<2e-2
        converged = true;
      elseif iter>2
        rmsdPctChange = abs(100*(rmsdNew-rmsdLast)/rmsdLast)
        converged = rmsdPctChange<50;
      end
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
    if ~isempty(Orients)
      gridPhi = repmat(Orients(1,:),[1,2^(iter-1)]);
      gridTheta = repmat(Orients(2,:),[1,2^(iter-1)]);
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

% average over trajectories for time-domain signal output
TDSignal = mean(TDSignal, 2);
TDSignal = TDSignal(1:Par.nSteps);

% Final processing
%-------------------------------------------------------------------------------

switch nargout
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
    title(sprintf('%0.8g GHz',Exp.mwFreq));
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

if iOrient==0, reverseStr = ''; return; end

avgTime = toc/iOrient;
secsLeft = (nOrient - iOrient)*avgTime;
minsLeft = floor(secsLeft/60);

secsElap = toc;
minsElap =  floor(secsElap/60);

msg3 = sprintf('  Time elapsed:   %02d:%02d:%02.0f (%g s/orientation)\n', floor(minsElap/60), mod(minsElap,60), mod(secsElap,60),avgTime);
msg4 = sprintf('  Time remaining: %02d:%02d:%02.0f\n', floor(minsLeft/60), mod(minsLeft,60), mod(secsLeft,60));
msg = [msg3, msg4];

if reverse
  fprintf([reverseStr, msg]);
else
  fprintf([msg]);  
end
reverseStr = repmat(sprintf('\b'), 1, length(msg));

end

function  y = zeropad(x, M)
  N = length(x);
  if iscolumn(x), y = [x; zeros(M-N, 1)]; end
  if isrow(x), y = [x, zeros(1, M-N)]; end
end
