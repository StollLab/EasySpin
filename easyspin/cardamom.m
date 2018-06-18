% cardamom  Trajectory-based simulation of CW-EPR spectra.
%
%   cardamom(Sys,Par,Exp)
%   cardamom(Sys,Par,Exp,Opt)
%   cardamom(Sys,Par,Exp,Opt,MD)
%   spc = cardamom(...)
%   [B,spc] = cardamom(...)
%   [B,spc,expectval,t] = cardamom(...)
%
%   Computes a CW-EPR spectrum of an 14N nitroxide radical using stochastic 
%   trajectories.
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
%     Coefs          numeric matrix, size = (nCoefs,2)
%                    array of orienting potential coefficients, with each row
%                    consisting of the corresponding real and imaginary parts
%
%     LMK            numeric matrix, size = (nCoefs,3)
%                    quantum numbers L, M, and K corresponding to each set of 
%                    coefficients
%
%     ProbDensFun    numeric array, 3D
%                    probability distribution grid to be used for
%                    calculating the pseudopotential and the torque
%
%     PseudoPotFun   numeric array, 3D
%                    orienting pseudopotential grid to be used for
%                    calculating the torque
%
%     TransRates     numeric matrix, size = (nStates,nStates)
%                    transition rate matrix describing inter-state dynamics
%                    for kinetic Monte Carlo simulations
%
%     TransProb      numeric matrix, size = (nStates,nStates)
%                    transition probability matrix describing inter-state 
%                    dynamics for kinetic Monte Carlo simulations
%
%     States         numeric matrix, size = (3,nStates)
%                    Euler angles for each state's orientation
%
%     Sys.lw         double or numeric vector, size = (1,2)
%                    vector with FWHM residual broadenings
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
%
% %     Sys.lwpp       double or numeric vector, size = (1,2)
% %                    peak-to-peak line widths, same format as Sys.lw
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
%     Omega          numeric, size = (3,1) or (3,nTraj)
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
%                    Brownian
%                    MOMD
%                    SRLS
%                    Discrete
%                    Molecular Dynamics
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
%                      the m_S=-1/2 subspace
%                    ISTOs: propagate the density matrix using
%                      irreducible spherical tensor operators
%
%    FFTWindow       1: use a Hamming window (default), 0: no window
%
%    truncate        double
%                    Time point (in nanoseconds) at which to stop using 
%                    full quantum dynamics propagator and begin using an
%                    approximate correlation function propagator.
%
%
%
%   MD: structure with molecular dynamics simulation parameters
%
%     dt             double
%                    time step for saving MD trajectory snapshots
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
%     GlobalDiff     double (optional)
%                    Diffusion coefficient for isotropic global rotational
%                    diffusion (s^-1)
%
%     TrajType       string
%                    Raw: raw MD trajectory data, including spin label, 
%                      protein, solvent, etc.
%                    Label: numeric array, trajectory of spin label atoms 
%                      only
%                    Frame: numeric array, trajectory of spin label 
%                      reference frame only
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
%                    Markov: coarse grain the trajectories by using the
%                      side chain dihedral angles to form a Markov state 
%                      model

%
%    Raw MD data input:
%
%     TrajFile       character array, or cell array containing character
%                    arrays as elements
%                    Name of trajectory output file(s), including the file
%                    extension ".[extension]".
%
%     AtomInfo
%
%         TopFile        character array
%                        Name of topology input file used for molecular 
%                        dynamics simulations.
%
%         ResName        character array
%                        Name of residue assigned to spin label side chain,
%                        e.g. "CYR1" is the default used by CHARMM-GUI.
%
%         AtomNames      structure array
%                        Structure array containing the atom names used in the 
%                        PSF to refer to the following atoms in the nitroxide 
%                        spin label molecule:
%
%                                   O (ONname)
%                                   |
%                                   N (NNname)
%                                  / \
%                        (C1name) C   C (C2name)
%
%   OR
%
%    Frame Trajectory input:
%
%     FrameX         numeric, size = (nSteps,3)
%                    x,y,z coordinate trajectory for X-axis vector
% 
%     FrameY         numeric, size = (nSteps,3)
%                    x,y,z coordinate trajectory for Y-axis vector
%  
%     FrameZ         numeric, size = (nSteps,3)
%                    x,y,z coordinate trajectory for Z-axis vector
%
%
%
%   Output:
%     B              numeric, size = (2*nSteps,1) 
%                    magnetic field (mT)
%
%     spc            numeric, size = (2*nSteps,1)
%                    derivative EPR spectrum
%
%     expectval      numeric, size = (2*nSteps,1)
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


function varargout = cardamom(Sys,Par,Exp,Opt,MD)

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
  case 3 % B,spc,expectval
  case 4 % B,spc,expectval,t
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

[Sys,err] = validatespinsys(Sys);
error(err);

if Sys.nElectrons>1, error('cardamom cannot be used for more than one electron.'); end
if Sys.nNuclei>1, error('cardamom cannot be used for more than one nucleus.'); end

if isfield(Sys, 'lw')
  if any(Sys.lw>0)
    isBroadening = 1;
  else
    isBroadening = 0;
  end
else
  isBroadening = 0;
end

% Check MD
% -------------------------------------------------------------------------

useMD = ~isempty(MD);

if useMD
  
  % scale the time axis
  if isfield(MD,'tScale')
    tScale = MD.tScale;
  else
    tScale = 1;
  end
  
  if ~isfield(MD,'dt')
    error('The MD trajectory time step MD.dt must be given.')
  end
  
  % check type of MD trajectory usage
  if ~isfield(MD,'TrajUsage')
    MD.TrajUsage = 'Explicit';
    MD.isFrame = 1;
  else
    if ~ischar(MD.TrajUsage)
      error('MD.TrajUsage must be a string.')
    end
    if ~any(strcmp({MD.TrajUsage},{'Explicit','Resampling','Markov'}))
      errmsg = sprintf('Entry ''%s'' for MD.TrajUsage not recognized.', MD.TrajUsage);
      error(errmsg)
    end
  end
  
  if ~isfield(MD,'TrajType')
    % use frame trajectories by default
    MD.TrajType = 'Frame';  % TODO add check for frame trajectory attributes
  else
    if ~any(strcmp({MD.TrajType},{'Frame','Dihedrals','Raw'}))
      errmsg = sprintf('Entry ''%s'' for MD.TrajType not recognized.', MD.TrajType);
      error(errmsg)
    end
  end
    
  switch MD.TrajType
    case 'Frame'
      if ~isfield(MD,'FrameX')||~isfield(MD,'FrameY')||~isfield(MD,'FrameZ')
        error('If using a frame trajectory, MD.FrameX, MD.FrameY, and MD.FrameZ must be given.')
      end
      sizeFrameX = size(MD.FrameX);

      if ~isequal(sizeFrameX,size(MD.FrameY))||~isequal(sizeFrameX,size(MD.FrameZ))
        error('All frame trajectory arrays in MD must have the same size.')
      end

      if sizeFrameX(2)~=3
        error('All frame trajectory arrays must be of size (nSteps,3).')
      end

      MD.nSteps = sizeFrameX(1);
    case 'Dihedrals'
      if ~strcmp(MD.TrajUsage,'Markov')
        error('Using a TrajType of Dihedrals requires TrajUsage to be Markov.')
      end
      
      if ~isfield(MD,'FrameX')||~isfield(MD,'FrameY')||~isfield(MD,'FrameZ')
        error('If using a frame trajectory, MD.FrameX, MD.FrameY, and MD.FrameZ must be given.')
      end
      sizeFrameX = size(MD.FrameX);

      if ~isequal(sizeFrameX,size(MD.FrameY))||~isequal(sizeFrameX,size(MD.FrameZ))
        error('All frame trajectory arrays in MD must have the same size.')
      end

      if sizeFrameX(2)~=3
        error('All frame trajectory arrays must be of size (nSteps,3).')
      end
      
      % use lag time to sample the trajectory such that that the result is
      % Markovian
      nLag = ceil(MD.tLag/MD.dt);
      
      dihedrals = [MD.chi1(1:nLag:end), ...
                   MD.chi2(1:nLag:end), ...
                   MD.chi4(1:nLag:end), ...
                   MD.chi5(1:nLag:end)];
                 
      Par.dt = tScale*MD.tLag;
    case 'Raw'
      if ~isfield(MD,'TrajFile')||~isfield(MD,'AtomInfo')
        error('MD.TrajFile and MD.AtomInfo must be provided.')
      end

      OutOpt.Verbosity = Opt.Verbosity;
      OutOpt.Format = 'Frame';  % TODO allow this to change depending on models chosen

      MD = mdload(MD.TrajFile, MD.AtomInfo, OutOpt);

      MD.dt = tScale*MD.dt;
      MD.nSteps = size(MD.FrameZ, 1);
    otherwise
      error('Entry for MD.TrajType not recognized.')
  end
  

  MD.FrameX = permute(MD.FrameX, [2, 3, 4, 1]);
  MD.FrameY = permute(MD.FrameY, [2, 3, 4, 1]);
  MD.FrameZ = permute(MD.FrameZ, [2, 3, 4, 1]);

  MDTrajLength = size(MD.FrameX, 4);

  MD.RTraj = zeros(3,3,1,MDTrajLength);
%   MD.RTraj(1,:,1,:) = MD.FrameX;
%   MD.RTraj(2,:,1,:) = MD.FrameY;
%   MD.RTraj(3,:,1,:) = MD.FrameZ;
  MD.RTraj(:,1,1,:) = MD.FrameX;
  MD.RTraj(:,2,1,:) = MD.FrameY;
  MD.RTraj(:,3,1,:) = MD.FrameZ;

%   q = rotmat2quat(MD.RTraj);
%   [alpha, beta, gamma] = quat2euler(q);

  % Check for orthogonality of rotation matrices
  RTrajInv = permute(MD.RTraj,[2,1,3,4]);

  if ~allclose(multimatmult(MD.RTraj,RTrajInv),...
               repmat(eye(3),1,1,size(MD.RTraj,3),size(MD.RTraj,4)),...
               1e-14)
    error('Rotation matrices obtained from frame trajectory are not orthogonal.')
  end

  MD.nTraj = size(MD.RTraj,3);

  
  if strcmp(MD.TrajUsage,'Markov')
    Opt.Model = 'Discrete';
    Opt.statesOnly = 1;
    
    MD.nStates = 48;
%       MD.nStates = 5;

    % Perform k-means clustering
    [MD.stateTraj,centroids] = clusterDihedrals(dihedrals,MD.nStates);
    MD.nSteps = size(MD.stateTraj, 1);  % TODO: find a way to process different step sizes here

    % Initialize HMM using clustering results
    mu0 = centroids.';
    for iState = 1:MD.nStates
      Sigma0(:,:,iState) = cov(dihedrals(MD.stateTraj==iState,:));
    end
%     prior0 = normalise(rand(MD.nStates,1));
%     transmat0 = mk_stochastic(rand(MD.nStates,MD.nStates));

    mixmat0 = ones(MD.nStates,1);

    randints = randi(size(MD.stateTraj,1), MD.nStates, 1);
    prior0 = MD.stateTraj(randints);
    transmat0 = calc_TPM(MD.stateTraj,MD.nStates).';
%     Sys.States0 =  TODO: find a way to assign the equilibrium distribution

    checkEmptyTrans = 1;
    ProbRatioThresh = 1e-3;  % threshold in probability ratio for finding 
                             % rarely visited states
    while checkEmptyTrans
      % run expectation-maximization algorithm on HMM model parameters
      ModelIn.prior = prior0;
      ModelIn.transmat = transmat0;
      ModelIn.mu = mu0;
      ModelIn.Sigma = Sigma0;

      [logL, ModelOut] = cardamom_emghmm(dihedrals, ModelIn, 20);

      prior1 = ModelOut.prior;
      transmat1 = ModelOut.transmat;
      mu1 = ModelOut.mu;
      Sigma1 = ModelOut.Sigma;
      
      nStates = size(transmat1,1);
      EmptyTransList = zeros(nStates,1);
      for iState = 1:nStates
        % check for rare states
        maxProb = max(transmat1(:));
        StateTransProbs = [transmat1(iState,:).'; transmat1(:,iState)];
        EmptyTransList(iState) = all(StateTransProbs./maxProb < ProbRatioThresh);
      end
      if any(EmptyTransList)
        % reset model vars by removing entries for rarely visited states
        idxNonEmpty = ~EmptyTransList;
        prior0 = prior1(idxNonEmpty);
        transmat0 = transmat1(idxNonEmpty,:);
        transmat0 = transmat0(:,idxNonEmpty);
        mu0 = mu1(:,idxNonEmpty);
        Sigma0 = Sigma1(:,:,idxNonEmpty);
        mixmat0 = mixmat1(idxNonEmpty);
      else
        % no rare states found, reset nStates if it has changed
        checkEmptyTrans = 0;
        MD.nStates = size(transmat1,1);
      end
    end
  end
  % estimate rotational diffusion time scale
  acorrX = autocorrfft(squeeze(MD.FrameX.^2), 2, 1, 1);
  acorrY = autocorrfft(squeeze(MD.FrameY.^2), 2, 1, 1);
  acorrZ = autocorrfft(squeeze(MD.FrameZ.^2), 2, 1, 1);

  N = round(MDTrajLength/4);

  % calculate correlation time
  time = linspace(0, N*MD.dt, N);
  tauRX = max(cumtrapz(time,acorrX(1:N)));
  tauRY = max(cumtrapz(time,acorrY(1:N)));
  tauRZ = max(cumtrapz(time,acorrZ(1:N)));
%   [k,c,yfit] = exponfit(time, acorr(1:N), 2, 'noconst');
%   tauR = 1/max(k);

%   DiffLocal = 1/6/(tauRZ);
  tauR = mean([tauRX, tauRY, tauRZ]);
  DiffLocal = 1/6/tauR;
  MD.tauR = tauR;

%   MD.FrameX = [];
%   MD.FrameY = [];
%   MD.FrameZ = [];
  RTrajInv = [];

end

% Check Par
% -------------------------------------------------------------------------

% If number of trajectories is not given, set it to 100
if ~isfield(Par,'nTraj')&&useMD==0, Par.nTraj = 100; end

% TODO add error checks from stochtraj and create a skipcheck flag for stochtraj

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
%   if (abs(Dt - tMax/nSteps) > eps), error('t is not linearly spaced.'); end
  
elseif isfield(Par,'nSteps') && isfield(Par,'dt')
    % number of steps and time step are given
    nStepsQuant = Par.nSteps;
    if ~isfield(Par,'Dt')
      Par.Dt = Par.dt;
    end
    if Par.Dt<Par.dt
      error('Par.dt must less than or equal to Par.Dt.')
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
  error(['You must specify a time array Par.t, or a number of steps '...
         'Par.nSteps and either a time step Par.dt or total time Par.tMax.'])
end

dtQuant = Par.Dt;
dtStoch = Par.dt;

% decide on a simulation model based on user input
if useMD
  if ~isfield(Par,'Model')
    % no Model given
    Par.Model = 'Molecular Dynamics';
  elseif ~strcmp(Par.Model,'Molecular Dynamics')
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
    Par.lag = ceil(3*MD.tauR/Par.Dt);  % use 2 ns lag between windows
    if Par.nSteps<Par.nBlocks
      % Par.nSteps not changed from user input
      Par.nTraj = floor((Par.nBlocks-Par.nSteps)/Par.lag) + 1;
    else
      Par.nSteps = nBlocks;
      Par.nTraj = 1;
    end
  end
  
else
  % no rotation matrices provided, so perform stochastic dynamics
  % simulations internally to produce them
  if ~isfield(Par,'Model')
    % no Model given
    if isfield(Sys,'LMK') && isfield(Sys,'Coefs') ...
       || isfield(Sys,'ProbDensFun') || isfield(Sys,'PseudoPotFun')
      % LMK and ordering coefs OR user-supplied potential given, so 
      % simulate MOMD
      Par.Model = 'MOMD';
    
    elseif xor(isfield(Sys,'LMK'),isfield(Sys,'Coefs'))
      error(['Both Sys.LMK and Sys.Coefs need to be declared for a MOMD '...
            'simulation.'])
    
    else
      % user did not specify a model or ordering potential, so perform 
      % Brownian simulation
      Par.Model = 'Brownian';
    end
  
  else
    % Model is specified
    if strcmp(Par.Model,'Brownian') && ...
        (isfield(Sys,'LMK')||isfield(Sys,'Coefs')||isfield(Sys,'ProbDensFun')||isfield(Sys,'PseudoPotFun'))
      error(['Conflicting inputs: Par.Model is set to "Brownian", but '...
            'Sys.LMK, Sys.Coefs, Sys.ProbDensFun, or Sys.PseudoPotFun have been declared.'])
    elseif strcmp(Par.Model,'MOMD') && (~isfield(Sys,'LMK')||~isfield(Sys,'Coefs'))
      if xor(isfield(Sys,'LMK'), isfield(Sys,'Coefs'))
        error('Both Sys.LMK and Sys.Coefs need to be declared for a MOMD simulation.')
      end
    elseif strcmp(Par.Model,'Molecular Dynamics') && (~isfield(MD,'RTraj')||~isfield(MD,'dt'))
      error('For Molecular Dynamics, both MD.RTraj and MD.dt need to be specified.')
    end
  end
  
  % check Par.Dt
  if Par.dt<Par.Dt
    Par.isBlock = 1;
  else 
    % same step size
    Par.isBlock = 0;
  end
    
  % find block properties for block averaging
  if Par.isBlock
    [Par.nBlocks,Par.BlockLength] = findblocks(Par.Dt, Par.dt, nStepsStoch);
  end
    
end

Model = Par.Model;

% Check Exp
% -------------------------------------------------------------------------

[Exp, CenterField] = validate_exp('cardamom',Sys,Exp);

omega = 2*pi*Exp.mwFreq*1e9;  % GHz -> rad s^-1

FieldSweep = true;  % TODO expand usage to include frequency sweep

% Set up horizontal sweep axis
xAxis = linspace(Exp.Range(1),Exp.Range(2),Exp.nPoints);  % field axis, mT


% Check Opt
% -------------------------------------------------------------------------

% if ~isfield(Opt,'chkcon'), chkcon = 0; end  % TODO implement spectrum convergence tests

if ~isfield(Opt,'Model')
  switch Par.Model
    case 'Brownian'
      Opt.Model = 'Continuous';
    case 'MOMD'
      Opt.Model = 'Continuous';
    case 'SRLS'
      Opt.Model = 'Continuous';
    case 'Discrete'
      Opt.Model = 'Discrete';
    case 'Molecular Dynamics'
      Opt.Model = 'Continuous';
  end
else
  switch Opt.Model
    case 'Continuous'
      if ~strcmp(Par.Model,'Brownian')||~strcmp(Par.Model,'MOMD')...
         ||~strcmp(Par.Model,'SRLS')||~strcmp(Par.Model,'Molecular Dynamics')
        error(['Continuous model is incompatible with the chosen Par.Model. '...
               'Please check documentation for Continuous model compatibilities.'])
      end
    case 'Discrete'
      if ~(strcmp(Par.Model,'Discrete')||strcmp(Par.Model,'Molecular Dynamics'))
        error(['Discrete model is incompatible with the chosen Par.Model. '...
               'Please check documentation for Discrete model compatibilities.'])
      end
  end
end

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

if ~isfield(Opt,'specCon')
  specCon = 0;
else
  specCon = Opt.specCon;
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

% Check dynamics and ordering
% -------------------------------------------------------------------------

Dynamics = validate_dynord('cardamom',Sys,FieldSweep);
Sys.Diff = Dynamics.Diff;

% Generate grids
% -------------------------------------------------------------------------

switch Model
  case 'Brownian'
    % no ordering present, so trajectory starting points are arbitrary
    if ~isfield(Par,'nOrients')
      % if Par.nOrients is not given, just use Par.nTraj as number of
      % orientations
      nOrients = Par.nTraj;
    else
      nOrients = Par.nOrients;
    end
    if specCon, nOrients = ceil(nOrients/2); end
    
  case 'MOMD'  %  TODO implement directors and ordering
    if ~isfield(Par,'nOrients')
      % if Par.nOrients is not given, just use Par.nTraj as number of
      % orientations
      nOrients = Par.nTraj;
    else
      nOrients = Par.nOrients;
    end
    
    if specCon
      nOrients = ceil(nOrients/2);
      skip = 0;
      gridPts = 2*sobol_generate(1,nOrients,skip)-1;
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    else
      if isfield(Par,'Orients')
        Orients = Par.Orients;
        gridPhi = Orients(1,:);
        gridTheta = Orients(2,:);
      else
        gridPts = linspace(-1,1,nOrients);
        gridPhi = sqrt(pi*nOrients)*asin(gridPts);
        gridTheta = acos(gridPts);
      end
    end
    
  case 'SRLS'  %  TODO implement multiple diffusion frames
    DiffLocal = Dynamics.Diff;
    if isfield(Sys, 'ProbDensFun')
      ProbDensFunLocal = Sys.ProbDensFun;
    end
    if isfield(Sys, 'PseudoPotFun')
      PseudoPotFunLocal = Sys.PseudoPotFun;
    end
    DiffGlobal = 6e6;
    if ~isfield(Par,'nOrients')
      % if Par.nOrients is not given, just use Par.nTraj as number of
      % orientations
      nOrients = Par.nTraj;
    else
      nOrients = Par.nOrients;
    end

    if specCon
      nOrients = ceil(nOrients/2);
      skip = 0;
      gridPts = 2*sobol_generate(1,nOrients,skip)-1;
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    else
      if isfield(Par,'Orients')
        Orients = Par.Orients;
        gridPhi = Orients(1,:);
        gridTheta = Orients(2,:);
      else
        gridPts = linspace(-1,1,nOrients);
        gridPhi = sqrt(pi*nOrients)*asin(gridPts);
        gridTheta = acos(gridPts);
      end
    end
    
  case 'Discrete'
    DiffLocal = Dynamics.Diff;
    DiffGlobal = 6e6;
    if ~isfield(Par,'nOrients')
      % if Par.nOrients is not given, just use Par.nTraj as number of
      % orientations
      nOrients = Par.nTraj;
    else
      nOrients = Par.nOrients;
    end

    if specCon
      nOrients = ceil(nOrients/2);
      skip = 0;
      gridPts = 2*sobol_generate(1,nOrients,skip)-1;
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    else
      if isfield(Par,'Orients')
        Orients = Par.Orients;
        gridPhi = Orients(1,:);
        gridTheta = Orients(2,:);
      else
        gridPts = linspace(-1,1,nOrients);
        gridPhi = sqrt(pi*nOrients)*asin(gridPts);
        gridTheta = acos(gridPts);
      end
    end
    
  case 'Molecular Dynamics' % TODO process RTraj based on size of input
    if ~isfield(Par,'nOrients')
      error('nOrients must be specified for the Molecular Dynamics model.')
    end
    nOrients = Par.nOrients;
    
    if specCon
      specCon, nOrients = ceil(nOrients/2);
      skip = 0;
      gridPts = 2*sobol_generate(1,nOrients,skip)-1;
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    else
      if isfield(Par,'Orients')
        Orients = Par.Orients;
        gridPhi = Orients(1,:);
        gridTheta = Orients(2,:);
      else
        gridPts = linspace(-1,1,nOrients);
        gridPhi = sqrt(pi*nOrients)*asin(gridPts);
        gridTheta = acos(gridPts);
      end
    end
    
    if strcmp(Opt.Method, 'ISTOs')
      % this method uses quaternions, not rotation matrices, so convert
      % MD.RTraj to quaternions here before the simulation loop
      MD.qTraj = rotmat2quat(MD.RTraj);
      clear MD.RTraj
    end

    if ~strcmp(MD.TrajUsage,'Explicit')
      switch MD.TrajUsage
        case 'Resampling'
          M = MDTrajLength;

          % calculate orienting potential energy function
          theta = squeeze(acos(MD.FrameZ(3,:,:,1:M)));
          phi = squeeze(atan2(MD.FrameY(3,:,:,1:M), MD.FrameX(3,:,:,1:M)));
          psi = squeeze(atan2(-MD.FrameZ(2,:,:,1:M), MD.FrameZ(1,:,:,1:M)));

%           % use eulang to obtain Euler angles
%           phi = zeros(MDTrajLength,1);
%           theta = zeros(MDTrajLength,1);
%           psi = zeros(MDTrajLength,1);
%           tic
%           for k = 1:MDTrajLength
%             [phi(k), theta(k), psi(k)] = eulang(MD.RTraj(:,:,1,k));
%             updateuser(k, MDTrajLength)
%           end

          nBins = 90;
          phiBins = linspace(-pi, pi, nBins);
          thetaBins = linspace(0, pi, nBins/2);
          psiBins = linspace(-pi, pi, nBins);

          [pdf, ~] = histcnd([phi,theta,psi], {phiBins,thetaBins,psiBins});
%            pdf = pdf/trapz(phiBins, trapz(thetaBins, trapz(psiBins, pdf)));

          pdf(end,:,:) = pdf(1,:,:);  % FIXME why does it truncate to zero in the phi direction?
          pdf = smoothn(pdf);
%           pdf = smooth3(pdf,'gaussian');
          save('pdf_smoothing.mat', 'pdf')
        case 'Markov'
          
      end
    end

  otherwise
    error('Model not recognized. Please check the documentation for acceptable models.')
    
end

logmsg(1,'-- time domain simulation -----------------------------------------');

logmsg(1, '-- Model: %s -----------------------------------------', Model);

logmsg(1, '-- Method: %s -----------------------------------------', Opt.Method);

% Run simulation
% -------------------------------------------------------------------------

clear cardamom_propagatedm

HistTot = 0;

converged = 0;
iOrient = 1;
iter = 1;
spcArray = [];
nOrientsTot = 0;

t = linspace(0, nStepsQuant*Par.Dt, nStepsQuant);

while ~converged
  tic
%   for iOrient = 1:nOrients
  % trajectories might differ in length, so we need cells for allocation
  ExpectVal = {};
  tCell = [];
  reverseStr = [];

  % temporary cells to store intermediate results
  iExpectVal = cell(1,nOrients);
  itCell = cell(1,nOrients);
%   while iOrient<nOrients+1
  for iOrient = 1:nOrients

  %   Par.Omega = [grid_phi(iOrient); grid_theta(iOrient); 0];

    % generate/process trajectories
    switch Model
  %     case 'Stochastic'
  %     case 'Molecular Dynamics'
      case 'Brownian'
        Sys.Diff = Dynamics.Diff;
        Par.dt = dtStoch;
        Par.nSteps = nStepsStoch;
        [trash, qTraj] = stochtraj(Sys,Par,Opt);
        Par.qTraj = qTraj;
        Par.RTraj = quat2rotmat(qTraj);
        Par.qLab = [];
        Par.RLab = [];

      case 'MOMD'
%         qMult = repmat(euler2quat(gridPhi(iOrient), gridTheta(iOrient), 0),...
%                        [1,Par.nTraj,Par.nSteps]);
        qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
                       [1,Par.nTraj,nStepsQuant]);
        
        Sys.Diff = Dynamics.Diff;
        Par.dt = dtStoch;
        Par.nSteps = nStepsStoch;
        [trash, qTraj] = stochtraj(Sys,Par,Opt);
        % generate quaternions for rotating to different grid points
        
        if strcmp(Opt.Method,'ISTOs')
          Par.qLab = qLab;
        else
          Par.RLab = quat2rotmat(qLab);
        end
        
        Par.qTraj = qTraj;
        Par.RTraj = quat2rotmat(qTraj);

      case 'SRLS'
        qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
               [1,Par.nTraj,nStepsQuant]);
        
        % local diffusion
        Sys.Diff = DiffLocal;
        if isfield(Sys, 'ProbDensFun')
          Sys.ProbDensFun = ProbDensFunLocal;
        end
        if isfield(Sys, 'PseudoPotFun')
          Sys.PseudoPotFun = PseudoPotFunLocal;
        end
        Par.dt = dtStoch;
        Par.nSteps = nStepsStoch;
        [trash, qTrajLocal] = stochtraj(Sys,Par,Opt);

        % global diffusion
        Sys.Diff = DiffGlobal;
        if isfield(Sys, 'ProbDensFun')
          Sys.ProbDensFun = [];
        end
        if isfield(Sys, 'PseudoPotFun')
          Sys.PseudoPotFun = [];
        end
        Par.dt = dtQuant;
        Par.nSteps = nStepsQuant;
        [trash, qTrajGlobal] = stochtraj(Sys,Par,Opt);
        qLab = quatmult(qLab,qTrajGlobal);
        
        Par.qTraj = qTrajLocal;
        Par.RTraj = quat2rotmat(qTrajLocal);
        Par.qLab = qLab;
        Par.RLab = quat2rotmat(qLab);
        
      case 'Discrete'
        qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
               [1,Par.nTraj,nStepsQuant]);
             
        % local diffusion
%         Sys.Diff = DiffLocal;
        Par.dt = dtStoch;
        Par.nSteps = nStepsStoch;
        Opt.Model = 'Discrete';
        [trash, qTrajLocal] = stochtraj(Sys,Par,Opt);
        
        % global diffusion
        Sys.Diff = DiffGlobal;
        Par.dt = dtQuant;
        Par.nSteps = nStepsQuant;
        Opt.Model = 'Continuous';
        [trash, qTrajGlobal] = stochtraj(Sys,Par,Opt);
        qLab = quatmult(qLab,qTrajGlobal);
        
        Par.qTraj = qTrajLocal;
        Par.RTraj = quat2rotmat(qTrajLocal);
        Par.qLab = qLab;
        Par.RLab = quat2rotmat(qLab);
        
      case 'Molecular Dynamics'
        
        switch MD.TrajUsage
          case 'Explicit'
            % powder grid rotation
  %           qMult = repmat(euler2quat(gridPhi(iOrient), gridTheta(iOrient), 0),...
  %                          [1,MD.nTraj,MD.nSteps]);
            qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
                           [1,Par.nTraj,nStepsQuant]);

            % global diffusion
            if isfield(MD, 'GlobalDiff')
              Sys.Diff = MD.GlobalDiff;
              Par.dt = dtQuant;
              Par.nSteps = nStepsQuant;  % TODO find a way to set this up with a separate time step properly
              [trash, qTrajGlobal] = stochtraj(Sys,Par,Opt);
              qLab = quatmult(qLab, qTrajGlobal);
            end

            if strcmp(Opt.Method,'ISTOs')
              Par.qLab = qLab;
            else
              Par.RLab = quat2rotmat(qLab);
            end

          case 'Resampling'
            qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
                           [1,Par.nTraj,Par.nSteps]);

            if ~isfield(Par,'Omega')
  %             % pick trajectory starting points by bootstrapping MD data
  %             randints = sort(randi(MD.nSteps, 1, Par.nTraj));
  %             Par.Omega = [phi(randints).'; theta(randints).'; psi(randints).'];
  %             [alphaSamples, betaSamples, gammaSamples] = rejectionsample3d(pdf, phiBins, thetaBins, psiBins, Par.nTraj);
  %             Par.Omega = [alphaSamples; 
  %                          betaSamples; 
  %                          gammaSamples];
            end

  %           [hist3D, dummy] = histcnd([alphaSamples.',betaSamples.',gammaSamples.'],...  
  %                                  {phiBins,thetaBins,psiBins});
  %                                
  %           save('hist3D.mat', 'hist3D')

            Sys.ProbDensFun = pdf;
            Sys.Diff = DiffLocal;
            Par.dt = dtStoch;
            Par.nSteps = nStepsStoch;
            [trash, qTraj] = stochtraj(Sys,Par,Opt);

            % global diffusion
            if isfield(MD, 'GlobalDiff')
              Sys.Diff = MD.GlobalDiff;
              Par.dt = dtQuant;
              Par.nSteps = nStepsQuant;
              [trash, qTrajGlobal] = stochtraj(Sys,Par,Opt);
              qLab = quatmult(qLab, qTrajGlobal);
            end

  %           qTraj = quatmult(qLab, qTraj);

            Par.qTraj = qTraj;
            Par.RTraj = quat2rotmat(qTraj);

            Par.qLab = qLab;
            Par.RLab = quat2rotmat(qLab);
            
          case 'Markov'
            qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
                           [1,Par.nTraj,nStepsQuant]);
%             qLab = repmat(euler2quat(0, gridTheta(iOrient), gridPhi(iOrient)),...
%                            [1,Par.nTraj]);

%             randints = sort(randi(size(MD.stateTraj,1), 1, Par.nTraj));
%             Sys.States0 = MD.stateTraj(randints).';  % FIXME sample from prior distribution
            Sys.States0 = rejectionsample(MD.nStates, prior1, Par.nTraj);
            Sys.TransProb = transmat1.';
            Par.dt = dtStoch;
            Par.nSteps = nStepsStoch;
            Opt.Model = 'Discrete';
            [trash, stateTraj] = stochtraj(Sys,Par,Opt);
            
           % global diffusion
            if isfield(MD, 'GlobalDiff')
              Sys.Diff = MD.GlobalDiff;
              Par.dt = dtQuant;
              Par.nSteps = nStepsQuant;
              Opt.Model = 'Continuous';
%               Par.Omega = qLab;
%               Par.Omega = [0;0;0];
              [trash, qTrajGlobal] = stochtraj(Sys,Par,Opt);
%               qLab = qTrajGlobal;
              qLab = quatmult(qLab, qTrajGlobal);
            end
            
            Par.RTraj = MD.RTraj(:,:,1,1:nLag:end);
            Par.stateTraj = stateTraj;
            Par.qLab = qLab;
            Par.RLab = quat2rotmat(qLab);
        end

    end

    % propagate the density matrix
    Par.nSteps = nStepsQuant;
    Par.Dt = dtQuant;
    Par.dt = dtStoch;
    Sprho = cardamom_propagatedm(Sys,Par,Opt,MD,omega,CenterField);

    % average over trajectories
    if strcmp(Opt.debug.EqProp,'time')
      Sprho = squeeze(mean(Sprho,3));
    end

    iExpectVal{1,iOrient} = 0;
    % calculate the expectation value of S_{+}
    for k = 1:size(Sprho,1)
      iExpectVal{1,iOrient} = iExpectVal{1,iOrient} + squeeze(Sprho(k,k,:));  % take traces TODO try to speed this up using equality tr(A*B)=sum(sum(A.*B))
    end
%     ExpectVal{1,iOrient} = squeeze(rho(1,1,:,:)+rho(2,2,:,:)+rho(3,3,:,:));  % take traces TODO try to speed this up using equality tr(A*B)=sum(sum(A.*B))

    if Opt.Verbosity
      updateuser(iOrient,nOrients)
    end

    itCell{1,iOrient} = t.';
    
    iOrient = iOrient + 1;

  end

  % Trajectory averaging and statistics
%   ExpectVal = [ExpectVal, cellfun(@(x) mean(x,1).', iExpectVal, 'UniformOutput', false)];

  % Store simulations at new starting orientations from each iteration
  ExpectVal = cat(2, ExpectVal, iExpectVal);
  tCell = cat(2, tCell, itCell);

%   if iter==1
%     expectval = cellfun(@(x) mean(x,1).', iexpectval, 'UniformOutput', false);
%     tcell = itcell;
%   else
%     expectval = [expectval, cellfun(@(x) mean(x,1).', iexpectval, 'UniformOutput', false)];
%     tcell = [tcell, itcell];
%   end
  
% end

%   if strcmp(Model, 'Molecular Dynamics')
%     % these variables can take up a lot of memory and might prevent the user 
%     % from implementing a fine enough grid for powder averaging 
%     clear RTraj
%     clear RTrajInv
%     clear qmult
%     clear MD.RTraj
%     clear Traj
%   end

  RTraj = [];
  qTraj = [];
  Par.RTraj = [];
  Par.qTraj = [];

% Perform FFT
% -------------------------------------------------------------------------

  % windowing
  if FFTWindow
  %   hamm = 0.54 + 0.46*cos(pi*t/max(t));
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
    M = ceil(tReq/Par.Dt);  % TODO make flexible for propagation length extensions
  else
  %   M = length(expectval);
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
      ExpectVal = exp(-alpha*tLong.^2).*ExpectVal;
  %     GaussBroad = cellfun(@(x) exp(-x.^2/TG.^2/8), tcell, 'UniformOutput', false);
  %     expectval = cellfun(@times, expectval, GaussBroad, 'UniformOutput', false);
    end
    if numel(Sys.lw)==2
      % Lorentzian broadening
      TL = Dynamics.T2; 
      ExpectVal = exp(-tLong/TL).*ExpectVal;
  %     LorentzBroad = cellfun(@(x) exp(-x/TL), tcell, 'UniformOutput', false);
  %     expectval = cellfun(@times, expectval, LorentzBroad, 'UniformOutput', false);
    end
  end

  % expectDt = cellfun(@times, expectval, tcell, 'UniformOutput', false);

  % spc = reshape(cell2mat(cellfun(@(x) fft(x,M), expectDt, 'UniformOutput', false)),...
  %               [M,nOrients]);
  
  % Multiply by t for differentiation and take the FFT
  spcArray = cat(2, spcArray, imag(fftshift(fft(ExpectVal.*tLong, [], 1))));
  spcNew = mean(spcArray,2);

  if specCon
    if iter==1
      spcLast = spcNew;
%       spread = std(spcArray,[],2);
      rmsdNew = 1;
    else
      span = max(spcNew)-min(spcNew);
      rmsdNew = sqrt(mean((spcNew-spcLast).^2))/span
      if rmsdNew<5e-4
        converged = 1;
      else
        rmsdPctChange = abs((rmsdNew-rmsdLast)/rmsdLast)
        converged = rmsdPctChange<10e-2;
      end
    end
%     if iter == 1
%       runningAvg = [];
%     end
%     gr = grstat(real(ExpectValArray));
%     converged = all(gr(:)<1.1);
  else
    converged = 1;
  end
  
  if converged
    spcAvg = spcNew;
    minsTot = floor(toc/60);
    msg = sprintf('Done!\nTotal simulation time: %d:%2.0f\n',minsTot,mod(toc,60));
    if Opt.Verbosity
      fprintf(msg);
    end
  else
    % increase the total number of orientations by 20%
    msg = sprintf('Convergence not achieved. Propagation is being extended.\n');
    if Opt.Verbosity
      fprintf(msg);
    end
    rmsdLast = rmsdNew;
    spcLast = spcNew;  % store for comparison after next iteration completes
    nOrientsTot = nOrientsTot + nOrients;
    skip = iter*nOrientsTot;  % seed Sobol sequence generator for next iteration
    
%     nOrients = ceil(0.2*nOrientsTot);  % simulate using 20% additional orientations
    nOrients = nOrientsTot;
    gridPts = 2*sobol_generate(1,nOrients,skip)-1;
    gridPhi = sqrt(pi*nOrients)*asin(gridPts);
    gridTheta = acos(gridPts);
    iter = iter + 1;
  end

end

freq = 1/(Par.Dt*M)*(-M/2:M/2-1);  % TODO check for consistency between FieldSweep and FFT window/resolution

% center the spectrum around the isotropic component of the g-tensor
fftAxis = mhz2mt(freq/1e6+Exp.mwFreq*1e3, mean(Sys.g));  % Note use of g0, not ge

% interpolate over horizontal sweep range
outspc = interp1(fftAxis, spcAvg, xAxis);

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

function [stateTraj,centroids] = clusterDihedrals(dihedrals,nStates)

chi1 = dihedrals(:,1);
chi2 = dihedrals(:,2);
% chi3 = dihedrals(:,3);
chi4 = dihedrals(:,3);
chi5 = dihedrals(:,4);

dihedrals = [wrapTo2Pi(chi1), wrapTo2Pi(chi2), wrapTo2Pi(chi4), chi5];

useParallel = false;
nReplicates = 5;
maxIter = 200;

% initialize cluster centroids
% chi1Min = wrapTo2Pi([-60;65;180]/180*pi);
% chi2Min = wrapTo2Pi([75;180]/180*pi);
% chi4Min = wrapTo2Pi([75;8;-100]/180*pi);
% chi5Min = wrapTo2Pi([180;77]/180*pi);

% start = zeros();

% opts = statset('Display', 'final', ...
%                'MaxIter', maxIter);
%                'UseParallel', useParallel);
% opts = statset('Display','final','MaxIter',maxIter,'UseParallel',useParallel,'Start',start);
% [stateTraj,centroids] = kmeans(dihedrals, nStates, ...
%                                'Distance', 'sqeuclidean', ...
%                                'Replicates', nReplicates, ...
%                                 'Options', opts);
[stateTraj,centroids] = cardamom_kmeans(dihedrals, nStates, 20, 1);

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

msg1 = sprintf('Iteration: %d/%d\n', iOrient, nOrient);
if avgTime<1.0
  msg2 = sprintf('%2.1f it/s\n', 1/avgTime);
else
  msg2 = sprintf('%2.1f s/it\n', avgTime);
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
    xSamples(iSample) = xGrid(xProposal);
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