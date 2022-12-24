% blochsteady   Steady-state solutions of Bloch equations
%
%  blochsteady(g,T1,T2,DeltaB0,B1,modAmp,modFreq)
%  blochsteady(g,T1,T2,DeltaB0,B1,modAmp,modFreq,Options)
%  [t,My] = blochsteady(...)
%  [t,Mx,My,Mz] = blochsteady(...)
%
%  Computes periodic steady-state solution of the Bloch equations
%  for a single spin-1/2 in the presence of a sinusoidal field modulation.
%
%  Inputs:
%    g        g value of the electron spin (S = 1/2)
%    T1       longitudinal relaxation time constant, µs
%    T2       transverse relaxation time constant, µs
%
%    DeltaB0  offset from resonance field, mT
%    B1       microwave field amplitude, mT
%    modAmp   peak-to-peak modulation field amplitude, mT
%    modFreq  modulation frequency, kHz
%
%    Options  calculation options
%      .Verbosity   whether to print information (0 or 1; 0 default)
%      .nPoints     number of points, chosen automatically by default
%      .kmax        highest Fourier order, chosen automatically by default
%      .Method      calculation method for time-domain signal
%                   'td'   explicit evolution in time-domain
%                   'fft'  using inverse Fourier transform (default)
%
%  Outputs:
%    t        time axis, µs
%    Mx       in-phase signal (dispersion)
%    My       quadrature signal (absorption)
%    Mz       longitudinal magnetization

% see
%   M. Tseitlin, G. R. Eaton, S. S. Eaton
%   Appl. Magn. Reson. 2013, 44, 1373-1379
%   https://doi.org/10.1007/s00723-013-0494-2

function varargout = blochsteady(g,T1,T2,DeltaB0,B1,modAmp,modFreq,Options)

if nargin==0, help(mfilename); return; end

switch nargin
  case 1, error('Second input (T1) is missing.');
  case 2, error('Third input (T2) is missing.');
  case 3, error('Fourth input (deltaB0) is missing.');
  case 4, error('Fifth input (B1) is missing.');
  case 5, error('Sixth input (ModAmp) is missing.');
  case 6, error('Seventh input (ModFreq) is missing.');
  case 7, Options = struct;
  case 8 % ok
  otherwise, error('Too many input arguments.');
end

if ~isstruct(Options)
  error('Eighth input (Options) must be a structure.');
end

if ~isfield(Options,'Verbosity'), Options.Verbosity = 0; end
global EasySpinLogLevel %#ok
EasySpinLogLevel = Options.Verbosity;

logmsg(1,['=begin=blochsteady======' char(datetime) '=================']);

if ~isfield(Options,'nPoints')
  Options.nPoints = [];
end
if ~isfield(Options,'Method')
  Options.Method = 'fft';
end

onlyAbsorption = nargout==2;

% Unit conversions
T1 = T1*1e-6; % µs -> s
T2 = T2*1e-6; % µs -> s

DeltaB0 = DeltaB0*1e-3; % mT -> T
B1 = B1*1e-3; % mT -> T
modAmp = modAmp*1e-3; % mT -> T
modFreq = modFreq*1e3; % kHz -> Hz
omegam = 2*pi*modFreq; % modulation angular frequency; Hz -> rad/s

M0 = 1;  % equilibrium magnetization
gamma = bmagn/hbar*g; % gyromagnetic ratio

% Some input checks
if numel(g)~=1
  error('g (1st input) must be a single number.');
end
if numel(T1)~=1 || T1<=0
  error('T1 (2nd input) must be a single positive number.');
end
if numel(T2)~=1 || T2<=0
  error('T2 (3rd input) must be a single positive number.');
end
if T2>T1
  error('T2 cannot be larger than T1.');
end
if numel(DeltaB0)~=1
  error('The field offset (4th input_ must be a single number.');
end
if numel(B1)~=1 || B1<=0
  error('The microwave field amplitude B1 (5th input) must be a single positive number.');
end
if numel(modAmp)~=1 || modAmp<=0
  error('The modulation amplitude (6th input) must be a single positive number.');
end
if numel(modFreq)~=1 || modFreq<=0
  error('The modulation frequency (7th input) must be a single positive number.');
end

logmsg(1,'Determination of maximum Fourier order (kmax)');
if ~isfield(Options,'kmax') || isempty(Options.kmax)
  % Estimator for maximum relevant Fourier order
  % (1) based on maximum field offset
  maxfieldoffset = max(modAmp/2-DeltaB0,DeltaB0+modAmp/2);
  maxfieldoffset = min(maxfieldoffset,modAmp);
  maxfreqoffset = (gamma/2/pi)*maxfieldoffset;
  kmax = maxfreqoffset/modFreq;
  kmax = ceil(kmax*1.4);
  logmsg(1,'  based on max. field offset: %d',kmax);
  
  % (2) based on T2
  %envelope = exp(-2/gam/ModAmp/T2*abs(fourierorder));
  threshold = 1e-6;
  thresholdorder = -log(threshold)/2*gamma*modAmp*T2;
  thresholdorder = ceil(thresholdorder);
  logmsg(1,'  based on T2: %d',thresholdorder);
  
  % Take the larger of (1) and (2), but at least some minimal k
  minkmax = 20;
  kmax = max(kmax,thresholdorder);
  kmax = max(kmax,minkmax);
  logmsg(1,'  larger of the two, but at least %d: %d',minkmax,kmax);
  
else
  kmax = Options.kmax; % largest Fourier order, |k|
  logmsg(1,'  user-supplied: %d',kmax);
end


% Solve Bloch equation for steady-state in frequency domain
%--------------------------------------------------------------------------
logmsg(1,'Frequency-domain steady-state solution');
logmsg(1,'  Set up diagonals');
k = (-(kmax+1):kmax+1).'; % use max. order kmax+1 to evaluate all terms
a = 1i*k*omegam;
b = gamma*DeltaB0;
c = gamma*B1;
d = gamma*modAmp/4;
tau1 = 1./(a+1/T1);
tau2 = 1./(a+1/T2);

q = (1:2*kmax+1)+1; % corresponds to -kmax:kmax range
c2m = d^2*tau2(q-1);
c2p = d^2*tau2(q+1);
c1m = b*d*(tau2(q) + tau2(q-1));
c1p = b*d*(tau2(q) + tau2(q+1));
c0 = a(q) + 1/T2 + b^2*tau2(q) + c^2*tau1(q) + d^2*(tau2(q-1) + tau2(q+1));
cL = c*M0*tau1(kmax+2)/T1;

% Assemble pentadiagonal coefficient matrix and RHS vector and solve sparse
% linear system of equations P*Y = C0 to get the Fourier coefficients Y = Mky
logmsg(1,'  Assemble sparse pentadiagonal matrix');
nRows = 2*kmax+1;
C0 = sparse(nRows,1);
C0(kmax+1) = cL;
P = spdiags([c2m c1m c0 c1p c2p],[0 1 2 3 4],nRows,nRows+4);
P = P(:,3:end-2);

logmsg(1,'  Solve sparse linear system for y Fourier coefficients');
logmsg(1,'    matrix size: %d x %d',size(P,1),size(P,2));
Y = P\C0;

% Calculate Mkx and Mkz from Mky
logmsg(1,'  Calculate x and z Fourier coefficients');
q = 2:2*kmax;  % drop one order (max order now is kmax-1)
if ~onlyAbsorption
  Xk = tau2(q+1).*(b*Y(q) + d*(Y(q-1) + Y(q+1)));
  deltak0 = zeros(2*kmax-1,1);
  deltak0(kmax) = 1;
  Zk = tau1(q+1).*(-c*Y(q) + M0*deltak0/T1);
end
Yk = Y(q);

% Sparse-to-full conversion (since ifft does not support sparse)
Yk = full(Yk);
if ~onlyAbsorption
  Xk = full(Xk);
  Zk = full(Zk);
end

% Compute time evolution of components
%--------------------------------------------------------------------------
logmsg(1,'Calculation of time-domain signal.');
tPeriod = 1/modFreq;  % modulation period

% Number of points in time domain
nPoints = Options.nPoints;
if isempty(nPoints)
  nPoints = 2*kmax-1;
end

t = linspace(0,tPeriod,nPoints+1).';
t(end) = [];

switch Options.Method
  case 'td'  % time-domain, explicit construction with complex exponentials
    logmsg(1,'  Method: explicit evolution with complex exponentials');
    if onlyAbsorption
      My = 0;
      for k = -(kmax-1):kmax-1
        idx = k+kmax; % index into array based on -(kmax-1):kmax-1
        phase = exp(1i*k*omegam*t);
        My = My + Yk(idx)*phase;
      end
      My = real(My); % remove numerical noise in imaginary part
    else
      Mx = 0;
      My = 0;
      Mz = 0;
      for k = -(kmax-1):kmax-1
        idx = k+kmax; % index into array based on -(kmax-1):kmax-1
        phase = exp(1i*k*omegam*t);
        Mx = Mx + Xk(idx)*phase;
        My = My + Yk(idx)*phase;
        Mz = Mz + Zk(idx)*phase;
      end
      f = @(x)max(abs(imag(x)))/max(abs(real(x)));
      f(Mx), f(My), f(Mz)
      Mx = real(Mx);
      My = real(My);
      Mz = real(Mz);
    end
    
  case 'fft'  % compute time domain signals using inverse Fourier transformation
    logmsg(1,'  Method: Inverse Fourier transform');
    logmsg(1,'    nPoints: %d, Mky: %d',nPoints,numel(Yk));
    n = numel(Yk);
    if nPoints<=n
      % Fourier transform
      My = n*real(ifft(ifftshift(Yk),'symmetric'));
      if ~onlyAbsorption
        Mx = n*real(ifft(ifftshift(Xk),'symmetric'));
        Mz = n*real(ifft(ifftshift(Zk),'symmetric'));
      end
      if nPoints<n
        tt = linspace(0,tPeriod,n+1);
        tt(end) = [];
        My = interp1(tt,My,t);
        if ~onlyAbsorption
          Mx = interp1(tt,Mx,t);
          Mz = interp1(tt,Mz,t);
        end
      end
    else
      fill = @(m) [m(1:kmax); zeros(nPoints-2*kmax+1,1); m(kmax+1:end)];
      FT = @(m)nPoints*real(ifft(fill(ifftshift(m)),'symmetric'));
      My = FT(Yk);
      if ~onlyAbsorption
        Mx = FT(Xk);
        Mz = FT(Zk);
      end
    end
    
  otherwise
    error('Unknown method.');
end

% Add last point if plotting
addLastPoint = nargout==0;
if addLastPoint
  t(end+1) = 1/modFreq;
  My(end+1) = My(1);
  if ~onlyAbsorption
    Mx(end+1) = Mx(1);
    Mz(end+1) = Mz(1);
  end
end

t = t/1e-6;    % seconds -> microseconds, for output and plotting

% Plotting
%--------------------------------------------------------------------------
switch nargout
  case 2
    varargout = {t,My};
  case 4
    varargout = {t,Mx,My,Mz};
  case 7
    varargout = {t,Mx,My,Mz,Xk,Yk,Zk};
  case 0
    varargout = {};
    plotresults;
end

logmsg(1,'=end=blochsteady========%s=================\n',char(datetime));

clear global EasySpinLogLevel

  function plotresults()

    xlineAvailable = ~verLessThan('matlab','9.5');  % R2018b

    tx1 = acos(-DeltaB0*2/modAmp)/omegam*1e6;
    tx2 = 1/modFreq*1e6 - tx1;
    reflineColor = [1 1 1]*0.7;
    
    clf
    subplot(3,1,[1,2]);
    plot(t,My);
    axis tight
    ylabel('absorption signal');
    title(sprintf('T_1 = %g µs, T_2 = %g µs, modAmp = %g mT, modFreq =  %g kHz, B_1 = %g mT, \\DeltaB = %g mT',...
      T1*1e6,T2*1e6,modAmp*1e3,modFreq/1e3,B1*1e3,DeltaB0*1e3));
    if xlineAvailable
      xline([tx1 tx2],'Color',reflineColor);
    end

    subplot(3,1,3);
    Bmod = modAmp/2*1e3*cos(omegam*t*1e-6);
    plot(t,Bmod);
    axis tight
    xlabel('time (µs)');
    ylabel('modulation field (mT)');
    if xlineAvailable
      yline(-DeltaB0*1e3,'Color',reflineColor,'Label','-\DeltaB0');
      xline([tx1 tx2],'Color',reflineColor);
    end
  end

end
