% cardamom  Trajectory-based simulation of CW-EPR spectra.
%
%   cardamom(Sys,Par,Exp)
%   cardamom(Sys,Par,Exp,Opt)
%   cardamom(Sys,Par,Exp,Opt,MD)
%   spc = cardamom(...)
%   [B,spc] = cardamom(...)
%   [B,spc,expectval] = cardamom(...)
%
%   Computes a CW-EPR spectrum of an 14N nitroxide radical using stochastic 
%   trajectories.
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
%     PseudoPotFun   numeric
%                    orienting pseudopotential function to be used for
%                    calculating the torque
%
%     Sys.lw         double or numeric, size = (1,2)
%                    vector with FWHM residual broadenings
%                         1 element:  GaussianFWHM
%                         2 elements: [GaussianFWHM LorentzianFWHM]
% %     Sys.lwpp       double or numeric, size = (1,2)
% %                    peak-to-peak line widths, same format as Sys.lw
%
%
%   Par: structure with simulation parameters
%     dt             double
%                    propagation time step (in seconds)
%
%     nSteps         int
%                    number of time steps per simulation
%
%     nTraj          int
%                    number of trajectories
%
%     Omega          numeric, size = (3,1) or (3,nTraj)
%                    Euler angles for starting orientation(s)
%
%     Model          string
%                    Brownian
%                    MOMD
%                    SRLS
%                    Molecular Dynamics
%
%
%   Exp: experimental parameter settings
%     mwFreq         double
%                    microwave frequency, in GHz (for field sweeps)
%
% %     Range          numeric, size = (1,2)
% %                    sweep range, [sweepmin sweepmax], in mT (for field sweep)
% %
% %     CenterSweep    numeric, size = (1,2)
% %                    sweep range, [center sweep], in mT (for field sweep)
% % 
% %     nPoints        int
% %                    number of points
% % 
% %     Harmonic       int
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
%                    Resampling: 
%
%    FFTWindow       1: use a Hamming window (default), 0: no window
%
%
%
%   MD: structure with molecular dynamics simulation parameters
%
%     tScale         double (optional
%                    scale the time step of the simulation based on
%                    incorrect diffusion coefficients, e.g. molecules 
%                    solvated in TIP3P water have diffusion coefficients
%                    that are ~2.5x too large
%
%     GlobalDiff     double (optional)
%                    Diffusion coefficient for isotropic global rotational
%                    diffusion (s^-1)
%
%     isFrame        integer
%                    0: raw MD data input
%                    1: frame trajectory input (MD data has already been 
%                       processed by mdload)
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
%                                   O (OName)
%                                   |
%                                   N (NName)
%                                  / \
%                        (C1Name) C   C (C2Name)
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
%                    expectation value of z-magnetization, 
%                    \langle S_{z} \rangle


function varargout = cardamom(Sys,Par,Exp,Opt,MD)
% Preprocessing
% -------------------------------------------------------------------------

switch nargin
  case 0
    help(mfilename); return;
  case 3 % Opt and MD not specified, initialize them
    Opt = struct;
    MD = struct;
  case 4 % MD not specified
    MD = struct;
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

if isfield(Sys, 'lw')
  if any(Sys.lw>0)
    Broadening = 1;
  else
    Broadening = 0;
  end
else
  Broadening = 0;
end

% Check MD
% -------------------------------------------------------------------------

useMD = ~isempty(fieldnames(MD));

if isfield(MD,'tScale')
  tScale = MD.tScale;  % diffusion constants of molecules solvated in TIP3P
                       % water are known to be too large
else
  tScale = 1;
end

if useMD
  if ~isfield(MD,'isFrame')
    % use frame trajectories by default
    MD.isFrame=1;  % TODO add check for frame trajectory attribute(s)
  end
  
  if MD.isFrame==1
    % use a frame trajectory
    if ~isfield(MD,'dt')
      error('The time step dt must be given in MD.')
    end
    
    if ~isfield(MD,'FrameX')||~isfield(MD,'FrameY')||~isfield(MD,'FrameZ')
      error('If using a frame trajectory, input MD must contain FrameX, FrameY, and FrameZ.')
    end
    sizeFrameX = size(MD.FrameX);

    if ~isequal(sizeFrameX,size(MD.FrameY))||~isequal(sizeFrameX,size(MD.FrameZ))
      error('All frame trajectory arrays in MD must have the same size.')
    end

    if sizeFrameX(2)~=3
      error('All frame trajectory arrays must be of size (nSteps,3).')
    end
    
    MD.nSteps = sizeFrameX(1);
  else
    % use raw MD input
    if ~isfield(MD,'TrajFile')||~isfield(MD,'AtomInfo')
%        ||~isfield(MD,'TopFile')...
%        ||~isfield(MD,'ResName')||~isfield(MD,'AtomNames')
      error('For MD input, TrajFile and AtomInfo must be provided.')
    end
    % generate rotation matrices from MD simulation data

    OutOpt.Verbosity = Opt.Verbosity;
    OutOpt.Frame = 1;

    MD = mdload(MD.TrajFile, MD.AtomInfo, OutOpt);
    
    MD.dt = tScale*MD.dt;
    MD.nSteps = size(MD.FrameZ, 1);
  end
  
  MD.FrameX = permute(MD.FrameX, [2, 3, 4, 1]);
  MD.FrameY = permute(MD.FrameY, [2, 3, 4, 1]);
  MD.FrameZ = permute(MD.FrameZ, [2, 3, 4, 1]);
  
  M = size(MD.FrameX, 4);
  
  MD.RTraj = zeros(3,3,1,M);
  MD.RTraj(:,1,1,:) = MD.FrameX;
  MD.RTraj(:,2,1,:) = MD.FrameY;
  MD.RTraj(:,3,1,:) = MD.FrameZ;
  
%   q = rotmat2quat(MD.RTraj);
%   [alpha, beta, gamma] = quat2euler(q);
  
  % Check for orthogonality of rotation matrices
  RTrajInv = permute(MD.RTraj,[2,1,3,4]);
  
  if ~allclose(matmult(MD.RTraj,RTrajInv),repmat(eye(3),1,1,size(MD.RTraj,3),size(MD.RTraj,4)),1e-14)
    error('The rotation matrices are not orthogonal.')
  end
  
  MD.nTraj = size(MD.RTraj,3);
  
  clear MD.FrameX
  clear MD.FrameY
  clear MD.FrameZ
end

% Check Par
% -------------------------------------------------------------------------

% If number of trajectories is not given, set it to 100
if ~isfield(Par,'nTraj')&&useMD==0, Par.nTraj = 100; end

% TODO add error checks from stochtraj and create a skipcheck flag for stochtraj

% decide on a simulation model based on user input
if useMD
  if ~isfield(Par,'Model')
    % no Model given
    Par.Model = 'Molecular Dynamics';
  elseif ~strcmp(Par.Model,'Molecular Dynamics')
    error('Mixing stochastic simulations with MD simulation input is not supported.')
  end
  
else
  % no rotation matrices provided, so perform stochastic dynamics
  % simulations internally to produce them
  if ~isfield(Par,'Model')
    % no Model given
    if isfield(Sys,'LMK') && isfield(Sys,'Coefs') ...
       || isfield(Sys, 'PseudoPotFun')
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
    if strcmp(Par.Model,'Brownian') && (isfield(Sys,'LMK')||isfield(Sys,'Coefs'))
      error(['Conflicting inputs: Par.Model is set to "Brownian", but at least '...
            'one of Sys.LMK or Sys.Coefs has been declared.'])
    elseif strcmp(Par.Model,'MOMD') && (~isfield(Sys,'LMK')||~isfield(Sys,'Coefs'))
      if xor(isfield(Sys,'LMK'), isfield(Sys,'Coefs'))
        error('Both Sys.LMK and Sys.Coefs need to be declared for a MOMD simulation.')
      end
    elseif strcmp(Par.Model,'Molecular Dynamics') && (~isfield(MD,'RTraj')||~isfield(MD,'dt'))
      error('For Molecular Dynamics, both MD.RTraj and MD.dt need to be specified.')
    end
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

if ~isfield(Opt,'Method')
  Opt.Method = 'Sezer';
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

% Check dynamics and ordering
% -------------------------------------------------------------------------

Dynamics = validate_dynord('cardamom',Sys,FieldSweep);

logmsg(1,'-- time domain simulation -----------------------------------------');


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
      gridPts = linspace(-1,1,nOrients);
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    end
    
  case 'SRLS'  %  TODO implement multiple diffusion frames
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
      gridPts = linspace(-1,1,nOrients);
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    end
    
  case 'Molecular Dynamics' % TODO process RTraj based on size of input
    if ~isfield(Par,'nOrients')
      error('nOrients must be specified for the Molecular Dynamics model.')
    end
    DiffGlobal = 6e6;
    nOrients = Par.nOrients;
    
    if specCon
      specCon, nOrients = ceil(nOrients/2);
      skip = 0;
      gridPts = 2*sobol_generate(1,nOrients,skip)-1;
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    else
      gridPts = linspace(-1,1,nOrients);
      gridPhi = sqrt(pi*nOrients)*asin(gridPts);
      gridTheta = acos(gridPts);
    end

    if strcmp(Opt.Method,'Resampling')
      
      % set up grid of starting orientations
%       if ~isfield(Par,'Omega')
%         Par.Omega = [sqrt(pi*Par.nTraj)*asin(linspace(0,1,Par.nTraj));...
%                      acos(linspace(0,1,Par.nTraj));...
%                      zeros(1,Par.nTraj)];
%       end
      
      % calculate orienting potential energy function
      theta = squeeze(acos(MD.FrameZ(3,:,:,:)));
      phi = squeeze(atan2(MD.FrameY(3,:,:,:), MD.FrameX(3,:,:,:)));
      psi = squeeze(atan2(-MD.FrameZ(2,:,:,:), MD.FrameZ(1,:,:,:)));

      nBins = 100;
      PhiBins = linspace(-pi, pi, nBins);
      ThetaBins = linspace(0, pi, nBins/2);
      PsiBins = linspace(-pi, pi, nBins);

      [PseudoPotFun, dummy] = histcnd([phi,theta,psi],...  
                                      {PhiBins,ThetaBins,PsiBins});
      
      PseudoPotFun(end,:,:) = PseudoPotFun(1,:,:);  % FIXME why does it truncate to zero in the phi direction?
      
      % first two dims are incorrectly ordered by histcnd
%       PseudoPotFun = permute(PseudoPotFun, [2, 1, 3]);  FIXME figure out correct dimension ordering

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% used for plotting PseudoPotFun on a sphere

% %       pad = (PseudoPotFun(1,:,:) + PseudoPotFun(end,:,:))/2;
%       PseudoPotFun(end,:,:) = PseudoPotFun(1,:,:);
%       
%       PseudoPotFun = smooth3(PseudoPotFun, 'gaussian');
% 
%       yy = permute(mean(PseudoPotFun, 3), [2, 1]);
%       yy = yy/max(yy(:));
%       yy(:,end) = yy(:,1);
%       
%       theta = linspace(0, pi, size(yy,1));                   % polar angle
%       phi = linspace(0, 2*pi, size(yy,2));                   % azimuth angle
%       
%       [Phi, Theta] = meshgrid(phi, theta);
%       radius = 1.0;
%       amplitude = 1.0;
% %       rho = radius + amplitude*yy;
%       rho = yy;
%       
%       r = radius.*sin(Theta);    % convert to Cartesian coordinates
%       x = r.*cos(Phi);
%       y = r.*sin(Phi);
%       z = radius.*cos(Theta);
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

      % estimate rotational diffusion time scale
      acorr = autocorrfft(squeeze(MD.FrameZ.^2), 1);

      N = round(length(MD.FrameZ)/2);
      M = round(N/2);

%       % the ACF will not always reach zero for a small number of trajectories,
%       % so subtract the offset, which does not change the correlation time
%       AutoCorrFFT = AutoCorrFFT - mean(AutoCorrFFT(M:3*M));

      % calculate correlation time
      time = linspace(0,N*MD.dt,N);
      tau = max(cumtrapz(time,acorr(1:N)));
%       tau = max(cumtrapz(time,datasmooth(acorr(1:N),500,'flat')));  % FIXME find a robust way to calculate this
      DiffLocal = 1/6/(2*tau);
%       Par.dt = tau/20;
    end

  otherwise
    error('Model not recognized. Please check the documentation for acceptable models.')
    
end

logmsg(1, '-- Model: %s -----------------------------------------', Model);

logmsg(1, '-- Method: %s -----------------------------------------', Opt.Method);

% Run simulation
% -------------------------------------------------------------------------

clear propagate_quantum

HistTot = 0;

converged = 0;
iOrient = 1;
iter = 1;
spcArray = [];
nOrientsTot = 0;

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
        
%         if strcmp(Opt.Method,'ISTOs')
%           % this method needs quaternions, not rotation matrices
%           [t, ~, qTraj] = stochtraj(Sys,Par,Opt);
%           Par.qTraj = qTraj;
%         else
%           % other methods use rotation matrices
%           [t, RTraj, ~] = stochtraj(Sys,Par,Opt);
%           Par.RTraj = RTraj;
%         end
        [t, RTraj, qTraj] = stochtraj(Sys,Par,Opt);
        Par.qTraj = qTraj;
        Par.RTraj = RTraj;

      case 'MOMD'
        [t, RTraj, qTraj] = stochtraj(Sys,Par,Opt);
        % generate quaternions for rotating to different grid points
        qMult = repmat(euler2quat(gridPhi(iOrient), gridTheta(iOrient), 0),...
                       [1,Par.nTraj,Par.nSteps]);
        qTraj = quatmult(qMult,qTraj);
        Par.qTraj = qTraj;
        Par.RTraj = quat2rotmat(qTraj);
%         if strcmp(Opt.Method,'ISTOs')
%           % this method needs quaternions, not rotation matrices
%           Par.qTraj = qTraj;
%         else
%           % other methods use rotation matrices
%           Par.RTraj = quat2rotmat(qTraj);
%         end

      case 'SRLS'
        Sys.Diff = DiffLocal;
        [t, ~, qTrajLocal] = stochtraj(Sys,Par,Opt);

        Sys.Diff = DiffGlobal;
        [t, ~, qTrajGlobal] = stochtraj(Sys,Par,Opt);
        qTraj = quatmult(qTrajGlobal,qTrajLocal);
        
        Par.qTraj = qTraj;
        Par.RTraj = quat2rotmat(qTraj);
%         if strcmp(Opt.Method,'ISTOs')
%           % this method needs quaternions, not rotation matrices
%           Par.qTraj = qTraj;  % ordering?
%         else
%           % other methods use rotation matrices
%           RTraj = quat2rotmat(qTraj);
%           Par.RTraj = RTraj;
%         end

      case 'Molecular Dynamics'
        % rotation matrices provided by external data, no need to do stochastic
        % simulation
        qMult = repmat(euler2quat(gridPhi(iOrient), gridTheta(iOrient), 0),...
                       [1,MD.nTraj,MD.nSteps]);
        MD.RTraj = matmult(quat2rotmat(qMult),MD.RTraj);

        if strcmp(Opt.Method,'Resampling')

          if ~isfield(Par,'Omega')
            % pick trajectory starting points by bootstrapping MD data
            randints = sort(randi(MD.nSteps, 1, Par.nTraj));
            Par.Omega = [phi(randints).'; theta(randints).'; psi(randints).'];
          end

          Sys.PseudoPotFun = PseudoPotFun;
          Sys.Diff = DiffLocal;
          [t, ~, qTraj] = stochtraj(Sys,Par,Opt);

          qMult = repmat(euler2quat(gridPhi(iOrient), gridTheta(iOrient), 0),...
                         [1,Par.nTraj,Par.nSteps]);
  %         qTraj = quatmult(qTrajGlobal, qTraj);
          qTraj = quatmult(qMult, qTraj);
  %         [alpha, beta, gamma] = quat2euler(qTraj);
  %         alpha = squeeze(alpha);
  %         beta = squeeze(beta);
  %         gamma = squeeze(gamma);
  %         
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
          
          Par.qTraj = qTraj;
          Par.RTraj = quat2rotmat(qTraj);
%           if strcmp(Opt.Method,'ISTOs')
%             % this method needs quaternions, not rotation matrices
%             Par.qTraj = qTraj;  % ordering?
%           else
%             % other methods use rotation matrices
%             Par.RTraj = quat2rotmat(qTraj);
%           end
        end

    end

    % propagate the density matrix
    rho = propagate_quantum(Sys,Par,Opt,MD,omega,CenterField);

    % average over trajectories
    rho = squeeze(mean(rho,3));

    % calculate the expectation value of S_{+}
    iExpectVal{1,iOrient} = squeeze(rho(1,1,:)+rho(2,2,:)+rho(3,3,:));  % take traces TODO try to speed this up using equality tr(A*B)=sum(sum(A.*B))
%     ExpectVal{1,iOrient} = squeeze(rho(1,1,:,:)+rho(2,2,:,:)+rho(3,3,:,:));  % take traces TODO try to speed this up using equality tr(A*B)=sum(sum(A.*B))

    if Opt.Verbosity
      updateuser(iOrient,nOrients)
    end

    if strcmp(Model,'Molecular Dynamics')
      nSteps = size(rho,4);
      t = linspace(0, nSteps*Par.dt, nSteps).';
    end

    itCell{1,iOrient} = t;
    
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
    hamm = cellfun(@(x) 0.54 + 0.46*cos(pi*x/max(x)), tCell, 'UniformOutput', false);
    ExpectVal = cellfun(@times, ExpectVal, hamm, 'UniformOutput', false);
  end

% zero padding for FFT to ensure sufficient B-field resolution (at most 0.1 G)
% expectval = cell2mat(cellfun(@(x) zeropad(x, maxlength), expectval, 'UniformOutput', false));
% tlong = (0:Par.dt:maxlength*Par.dt);

  Bres = 0.1; % G
  tReq = 1/(mt2mhz(Bres/10)*1e6); % mT -> s

  % if max(t)<treq
  tMax = max(cellfun(@(x) max(x), tCell));
  if tMax<tReq
    M = ceil(tReq/Par.dt);  % TODO make flexible for propagation length extensions
  else
  %   M = length(expectval);
    M = ceil(tMax/Par.dt);
  end
  ExpectVal = cell2mat(cellfun(@(x) zeropad(x, M), ExpectVal, 'UniformOutput', false));
  tLong = linspace(0, M*Par.dt, M).';

  % broadening
  if Broadening
    if Sys.lw(1)>0
      % Gaussian broadening
      TG = 1/(mt2mhz(Sys.lw(1))*1e6);
      ExpectVal = exp(-tLong.^2/TG^2/8).*ExpectVal;
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

freq = 1/(Par.dt*M)*(-M/2:M/2-1);  % TODO check for consistency between FieldSweep and FFT window/resolution

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
end

clear global EasySpinLogLevel
    
end

% Helper functions
% -------------------------------------------------------------------------

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

function  y = zeropad(x, M)
  N = length(x);
  if iscolumn(x), y = [x; zeros(M-N, 1)]; end
  if isrow(x), y = [x, zeros(1, M-N)]; end
end