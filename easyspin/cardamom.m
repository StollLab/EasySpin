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
% %     Sys.lw         double or numeric, size = (1,2)
% %                    vector with FWHM residual broadenings
% %                         1 element:  GaussianFWHM
% %                         2 elements: [GaussianFWHM LorentzianFWHM]
% %     Sys.lwpp       double or numeric, size = (1,2)
% %                    peak-to-peak line widths, same format as Sys.lw
%
%
%   Par: structure with simulation parameters
%     dt             double
%                    time step (in seconds)
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
%    Method           string
%                     'Sezer': propagate the density matrix using the 
%                      m_s=-1/2 subspace
%                     'DeSensi': propagate the density matrix using an 
%                      eigenvalue method
%
%
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
%   Sezer, et al., J.Chem.Phys. 128, 165106 (2008)
%     http://dx.doi.org/10.1063/1.2908075

function varargout = cardamom(Sys,Par,Exp,Opt)
% Preprocessing
% -------------------------------------------------------------------------

switch nargin
  case 0
    help(mfilename); return;
  case 3 % Opt not specified, initialize it
    Opt = struct;
  case 4 % Sys, Par, Exp, Opt specified
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

% decide on a simulation model based on user input
% if ~isfield(Sys,'LMK') && ~isfield(Sys,'Coefs')
  % no orienting potential specified, Brownian motion will be simulated
%   model = 'Brownian';
% else
  % 
%   model = 'MOMD';
% else
  %
%   model = 'SRLS';
% else
%   error('Input could not be interpreted based on the available models. Please check documentation.')
% end


% Check Par
% -------------------------------------------------------------------------

% If number of trajectories is not given, set it to 100
if ~isfield(Par, 'nTraj'), Par.nTraj = 100; end

dt = Par.dt;


% Check Exp
% -------------------------------------------------------------------------


[Exp, CenterField] = check_exp('cardamom',Sys,Exp);

Sys.B = CenterField/1000;  % TODO replace with omega0 or use fieldsweep
% Sys.B = mhz2mt(Exp.mwFreq*1e3, sum(Sys.g)/3)/1000


% Check Opt
% -------------------------------------------------------------------------

if ~isfield(Opt,'chkcon'), chkcon = 0; end


% Check dynamics and ordering
% -------------------------------------------------------------------------

FieldSweep = true;  % TODO fix this
Dynamics = check_dynord('cardamom',Sys,FieldSweep);

logmsg(1,'-- time domain simulation -----------------------------------------');


% Pre-allocate cells
% -------------------------------------------------------------------------

% we need cells here because the size of the dimension nSteps might differ
% between orientations due to differing times before convergence
if ~isfield(Par,'Model')
  Model='Brownian'; 
else
  Model=Par.Model;
end

switch Model
  case 'Brownian'
    logmsg(1,'-- model: Brownian dynamics -----------------------------------------');
    % no ordering present, so trajectory starting points are arbitrary
    if ~isfield(Par,'nOrients')
      nOrients = Par.nTraj;
    else
      nOrients = Par.nOrients;
    end
  case 'MOMD'  %  TODO implement directors and ordering
    logmsg(1,'-- model: MOMD -----------------------------------------');
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
end

expval = cell(1,nOrients);
tcell = cell(1,nOrients);


% Simulation
% -------------------------------------------------------------------------

tic
for iOrient = 1:nOrients

% Par.Omega = [grid_phi(iOrient); grid_theta(iOrient)];

[t, RTraj, qTraj] = stochtraj(Sys,Par);
tcell{1,iOrient} = t;
if strcmp(Model,'MOMD')
  qmult = repmat(euler2quat(0, grid_theta(iOrient), grid_phi(iOrient)),1,Par.nTraj,Par.nSteps);
  RTraj = quat2rotmat(quatmult(qmult,qTraj));
elseif strcmp(Model,'SRLS')
  
end

rho_t = propagate_quantum(Sys,Par,Opt,RTraj);
temp = squeeze(sum(rho_t,3));  % average over trajectories
expval{1,iOrient} = squeeze(temp(1,1,:)+temp(2,2,:)+temp(3,3,:));  % take traces TODO try to speed this up using equality tr(A*B)=sum(sum(A.*B))

updateuser(iOrient,nOrients)

end
mins_tot = floor(toc/60);
msg = sprintf('Done!\nTotal simulation time: %d:%2.0f\n',mins_tot,mod(toc,60));
fprintf(msg);
    

% Perform FFT

% hamm = 0.54 + 0.46*cos(pi*t/max(t));
% TL = 110e-9;  TODO implement Lorentzian and Gaussian broadening
TL = Dynamics.T2;

% Convolve with Lorentzian and multiply by t for differentiation
tdiff = cellfun(@(x) x.*exp(-x/TL), tcell, 'UniformOutput', false);
expvalDt = cellfun(@times, expval, tdiff, 'UniformOutput', false);

expval = mean(cell2mat(expval),2);

M = ceil(1e-6/dt);
spc = cell2mat(cellfun(@(x) fft(x,M), expvalDt, 'UniformOutput', false));
outspec = imag(fftshift(sum(spc,2)));
freq = 1/(dt*M)*(-M/2:M/2-1);  % TODO check for consistency between FieldSweep and FFT window/resolution

xAxis = freq;


% Final processing
% -------------------------------------------------------------------------

% borrowed from chili
switch (nargout)
case 0
  cla
%   plot(mhz2mt(xnum)+340, ynum)
  plot(xAxis,outspec)
  ylabel('Im(FFT(M_{+}(t)))')
  xlabel('f (MHz)')
%   if FieldSweep  TODO fix output plotting
%     if (xAxis(end)<10000)
%       plot(xAxis,outspec);
%       xlabel('magnetic field (mT)');
%     else
%       plot(xAxis/1e3,outspec);
%       xlabel('magnetic field (T)');
%     end
%     axis tight
%     ylabel('intensity (arb.u.)');
%     title(sprintf('%0.8g GHz, %d points',Exp.mwFreq,numel(xAxis)));
%   else
%     if (xAxis(end)<1)
%       plot(xAxis*1e3,spc);
%       xlabel('frequency (MHz)');
%     else
%       plot(xAxis,spc);
%       xlabel('frequency (GHz)');
%     end
%     axis tight
%     ylabel('intensity (arb.u.)');
%     title(sprintf('%0.8g mT, %d points',Exp.Field,numel(xAxis)));
%   end
case 1
  varargout = {outspec};
case 2
  varargout = {xAxis,outspec};
case 3
  varargout = {xAxis,outspec,expval};
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

if iOrient==nOrient
  clear reverseStr
end

end