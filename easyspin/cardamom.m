% cardamom  Trajectory-based simulation of CW-EPR spectra.
%
%   cardamom(Sys,Par,Exp,Opt)
%   spc = cardamom(...)
%   [B,spc] = cardamom(...)
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
%
% %     Sys.lw         double or numeric, size = (1,2)
% %                    vector with FWHM residual broadenings
% %                         1 element:  GaussianFWHM
% %                         2 elements: [GaussianFWHM LorentzianFWHM]
% %     Sys.lwpp       double or numeric, size = (1,2)
% %                    peak-to-peak line widths, same format as Sys.lw
%
%
%   MD: structure with molecular dynamics simulation parameters
%     RTraj          numeric, size = (3,3,N,M)
%                    rotation matrices, calculated externally, e.g. by
%                    performing a molecular dynamics simulation, where N is
%                    the number of trajectories and M is the number of time
%                    steps
%
%     dt             double
%                    MD simulation time step (in seconds)
%
% %     dcd           dcd file output from a molecular dynamics simulation
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
%     alpha          double or numeric, size = (1,nTraj)
%                    Euler angle alpha for starting orientation(s)
%
%     beta           double or numeric, size = (1,nTraj)
%                    Euler angle beta for starting orientation(s)
%
%     gamma          double or numeric, size = (1,nTraj)
%                    Euler angle gamma for starting orientation(s)
%                    
%
%
%
% %   Exp: experimental parameter settings
% %     mwFreq         double
% %                    microwave frequency, in GHz (for field sweeps)
% %
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
%     chkcon         if equal to 1, after the first nSteps of the 
%                    trajectories are calculated, both inter- and intra-
%                    trajectory convergence is checked using the Gelman-
%                    Rubin R statistic such that R<1.1, and if this 
%                    condition is not satisfied, then propagation will be 
%                    extended by either a length of time equal to the 
%                    average of tcorr or by 20% more time steps, whichever 
%                    is larger
%
%    Verbosity       0: no display, 1: show info
%
%    Method          string
%                    'Sezer': propagate the density matrix using an 
%                    analytical expression for the matrix exponential in 
%                    the m_s=-1/2
%                    'DeSensi': propagate the density matrix using an 
%                    eigenvalue method
%                    'Oganesyan': propagate the density matrix using
%                    irreducible spherical tensor operators and correlation 
%                    functions
%
%    FFTWindow       1: use a Hamming window (default), 0: no window
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
%     expval         numeric, size = (2*nSteps,1)
%                    expectation value of z-magnetization, 
%                    \langle S_{z} \rangle

% Implementation based on 
%   Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%     http://dx.doi.org/10.1063/1.2908075

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
  case 3 % B,spc,expval
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
EasySpinLogLevel = Opt.Verbosity;


% Check Sys
% -------------------------------------------------------------------------

[Sys,err] = validatespinsys(Sys);
error(err);

% Check MD
% -------------------------------------------------------------------------

useMD = ~isempty(fieldnames(MD));

if useMD
  if ~isfield(MD,'RTraj')||~isfield(MD,'dt')
    error('For Molecular Dynamics, both MD.RTraj and MD.dt need to be specified.')
  end
  % rotation matrices provided externally
  RTraj = MD.RTraj;
  if ~isnumeric(RTraj)||size(RTraj,1)~=3||size(RTraj,2)~=3||ndims(RTraj)~=4
    error('RTraj must be a 4D array of rotation matrices of size (3,3,MD.nTraj,MD.nSteps).')
  end
  % Check for orthogonality of rotation matrices
  RTrajInv = permute(RTraj,[2,1,3,4]);

  if ~allclose(matmult(RTraj,RTrajInv),repmat(eye(3),1,1,size(RTraj,3),size(RTraj,4)),1e-14)
    error('The rotation matrices are not orthogonal.')
  end
  
  MD.nTraj = size(RTraj,3);
  MD.nSteps = size(RTraj,4);
  
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
    error('Mixing stochastic simulations with external rotation matrices is not supported.')
  end
else
  % no rotation matrices provided, so perform stochastic dynamics
  % simulations internally to produce them
  if ~isfield(Par,'Model')
    % no Model given
    if isfield(Sys,'LMK') && isfield(Sys,'Coefs')
      % LMK and ordering coefs given, so simulate MOMD
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
      error('Both Sys.LMK and Sys.Coefs need to be declared for a MOMD simulation.')
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
  fftWindow = Opt.FFTWindow;
else
  fftWindow = 1;
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
  case 'MOMD'  %  TODO implement directors and ordering
    logmsg(1,'-- Model: MOMD -----------------------------------------');
    if ~isfield(Par,'nOrients')
      error('nOrients must be specified for the MOMD model.')
    end
    nOrients = Par.nOrients;
    grid_pts = linspace(-1,1,nOrients);
    grid_phi = sqrt(pi*nOrients)*asin(grid_pts);
    grid_theta = acos(grid_pts);
%   case 'SRLS'  %  TODO implement multiple diffusion frames
%     logmsg(1,'-- model: SRLS -----------------------------------------');  
%     
  case 'Molecular Dynamics' % TODO process RTraj based on size of input
    if ~isfield(Par,'nOrients')
      error('nOrients must be specified for the Molecular Dynamics model.')
    end
    nOrients = Par.nOrients;
    grid_pts = linspace(-1,1,nOrients);
    grid_phi = sqrt(pi*nOrients)*asin(grid_pts);
    grid_theta = acos(grid_pts);
  otherwise
    error('Model not recognized. Please check the documentation for acceptable models.')
end

logmsg(1, '-- Model: %s -----------------------------------------', Model);

logmsg(1, '-- Method: %s -----------------------------------------', Opt.Method);

% trajectories might differ in length, so we need cells for allocation
expect = cell(1,nOrients);
tcell = cell(1,nOrients);


% Run simulation
% -------------------------------------------------------------------------

clear updateuser
clear propagate_quantum

tic
for iOrient = 1:nOrients

  % Par.Omega = [grid_phi(iOrient); grid_theta(iOrient)];

  % generate/process trajectories
  switch Model 
    case 'Brownian'
      [t, RTraj, qTraj] = stochtraj(Sys,Par);
      if strcmp(Opt.Method,'Oganesyan')
        % this method needs quaternions, not rotation matrices
        Par.qTraj = qTraj;
      else
        % other methods use rotation matrices
        Par.RTraj = RTraj;
      end
    case 'MOMD'
      [t, RTraj, qTraj] = stochtraj(Sys,Par);
      qmult = repmat(euler2quat(grid_phi(iOrient), grid_theta(iOrient), 0),...
                     [1,Par.nTraj,Par.nSteps]);
      if strcmp(Opt.Method,'Oganesyan')
        % this method needs quaternions, not rotation matrices
        Par.qTraj = quatmult(qmult,qTraj);
      else
        % other methods use rotation matrices
        Par.RTraj = quat2rotmat(quatmult(qmult,qTraj));
      end
%   case 'SRLS'
    case 'Molecular Dynamics'
      % rotation matrices provided by external data, no need to do stochastic
      % simulation
      qmult = repmat(euler2quat(grid_phi(iOrient), grid_theta(iOrient), 0),...
                     [1,MD.nTraj,MD.nSteps]);
      MD.RTraj = matmult(quat2rotmat(qmult),RTraj);

  end

  % propagate the density matrix
  rho_t = propagate_quantum(Sys,Par,Opt,MD,omega,CenterField);
  
  % average over trajectories
  rho_t = squeeze(mean(rho_t,3));
  
  % calculate the expectation value of S_{+}
  expect{1,iOrient} = squeeze(rho_t(1,1,:)+rho_t(2,2,:)+rho_t(3,3,:));  % take traces TODO try to speed this up using equality tr(A*B)=sum(sum(A.*B))

  if Opt.Verbosity
    updateuser(iOrient,nOrients)
  end
  
  if strcmp(Model,'Molecular Dynamics')
    nSteps = size(rho_t,3);
    t = linspace(0, nSteps*Par.dt, nSteps).';
  end
  
  tcell{1,iOrient} = t;

end

mins_tot = floor(toc/60);
msg = sprintf('Done!\nTotal simulation time: %d:%2.0f\n',mins_tot,mod(toc,60));
if Opt.Verbosity
  fprintf(msg);
end

if strcmp(Model, 'Molecular Dynamics')
  % these variables can take up a lot of memory and might prevent the user 
  % from implementing a fine enough grid for powder averaging 
  clear RTraj
  clear RTrajInv
  clear qmult
  clear MD.RTraj
end

% Perform FFT
% -------------------------------------------------------------------------

if fftWindow
  hamm = 0.54 + 0.46*cos(pi*t/max(t));
else
  hamm = 1;
end
TL = Dynamics.T2;  % TODO implement Lorentzian and Gaussian broadening

% Convolve with Lorentzian and multiply by t for differentiation
tdiff = cellfun(@(x) hamm.*x.*exp(-x/TL), tcell, 'UniformOutput', false);
expectDt = cellfun(@times, expect, tdiff, 'UniformOutput', false);

% zero padding for FFT to ensure sufficient B-field resolution 
% (at most 0.1 mT)
Bres = 0.1; % mT
treq = 1/(mt2mhz(Bres)*1e6); % mT -> s
if max(t)<treq
  M = ceil(treq/Par.dt);  % TODO make flexible for propagation length extensions
else
  M = length(expect);
end

% relaxation and differentiation of spectrum via convolution
spc = reshape(cell2mat(cellfun(@(x) fft(x,M), expectDt, 'UniformOutput', false)),...
              [M,nOrients]);
spc = imag(fftshift(mean(spc,2)));
freq = 1/(Par.dt*M)*(-M/2:M/2-1);  % TODO check for consistency between FieldSweep and FFT window/resolution

% center the spectrum around the isotropic component of the g-tensor
fftAxis = mhz2mt(freq/1e6+Exp.mwFreq*1e3,mean(Sys.g));  % Note use of g0, not ge

% interpolate over horizontal sweep range
outspec = interp1(fftAxis,spc,xAxis);

% average over trajectories for <S_+(t)> output
tmat = cell2mat(tcell);
expect = mean(cell2mat(expect).*exp(-tmat/TL),2);

% Final processing
% -------------------------------------------------------------------------

switch (nargout)
case 0
  cla
  if FieldSweep  % TODO fix output plotting
    if (xAxis(end)<10000)
      plot(xAxis,outspec);
      xlabel('magnetic field (mT)');
    else
      plot(xAxis/1e3,outspec);
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
  varargout = {outspec};
case 2
  varargout = {xAxis,outspec};
case 3
  varargout = {xAxis,outspec,expect};
end

clear global EasySpinLogLevel
    
end

% Helper functions
% -------------------------------------------------------------------------

function updateuser(iOrient,nOrient)
% Update user on progress

persistent reverseStr

if isempty(reverseStr), reverseStr = ''; end

avg_time = toc/iOrient;
secs_left = (nOrient - iOrient)*avg_time;
mins_left = floor(secs_left/60);

msg1 = sprintf('Iteration: %d/%d\n', iOrient, nOrient);
if avg_time<1.0
  msg2 = sprintf('%2.1f it/s\n', 1/avg_time);
else
  msg2 = sprintf('%2.1f s/it\n', avg_time);
end
msg3 = sprintf('Time left: %d:%2.0f\n', mins_left, mod(secs_left,60));
msg = [msg1, msg2, msg3];

fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));

end