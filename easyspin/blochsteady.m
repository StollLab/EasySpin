% blochsteady   Steady-state solutions of Bloch equations
%
%  blochsteady(g,T1,T2,deltaB0,B1,ModAmp,ModFreq)
%  blochsteady(g,T1,T2,deltaB0,B1,ModAmp,ModFreq,Options)
%  [t,My] = blochsteady(...)
%  [t,Mx,My,Mz] = blochsteady(...)
%
%  Computes periodic steady-state solution of the Bloch equations
%  for a single spin-1/2 in the presence of field modulation.
%
%  Inputs:
%    g        g value of the electron spin (S = 1/2)
%    T1       longitudinal relaxation time, us
%    T2       transverse relaxation time, us
%
%    deltaB0  offset from resonance field, mT
%    B1       microwave field amplitude, mT
%    ModAmp   peak-to-peak field modulation amplitude, mT
%    ModFreq  modulation frequency, kHz
%
%  Outputs:
%    t        time axis, us
%    Mx       in-phase signal (dispersion)
%    My       quadrature signal (absorption)
%    Mz       longitudinal magnetization

function varargout = blochsteady(g,T1,T2,deltaB0,B1,ModAmp,ModFreq,Options)

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
global EasySpinLogLevel
EasySpinLogLevel = Options.Verbosity;


logmsg(1,['=begin=blochsteady======' datestr(now) '=================']);

if ~isfield(Options,'nPoints')
  Options.nPoints = 1000;
end
if ~isfield(Options,'Method')
  Options.Method = 'fft';
end

onlyAbsorption = nargout==2;

% Unit conversions
T1 = T1*1e-6; % us -> s
T2 = T2*1e-6; % us -> s

deltaB0 = deltaB0*1e-3; % mT -> T
B1 = B1*1e-3; % mT -> T
ModAmp = ModAmp*1e-3; % mT -> T
ModFreq = ModFreq*1e3; % kHz -> Hz

% Some range checks
if numel(T1)~=1 || (T1<=0)
  error('T2 must be a single nonnegative number.');
end
if numel(T2)~=1 || (T2<=0)
  error('T2 must be a single nonnegative number.');
end
if (T2>T1)
  error('T2 cannot be larger than T1.');
end
if numel(B1)~=1 || (B1<=0)
  error('The microwave field amplitude B1 must be a single nonnegative number.');
end
if numel(ModFreq)~=1 || (ModFreq<=0)
  error('The modulation frequency ModFreq must be a single nonnegative number.');
end
if numel(ModAmp)~=1 || (ModAmp<=0)
  error('The modulation amplitude ModAmp must be a single nonnegative number.');
end
if numel(deltaB0)~=1
  error('The field offset must be a single number.');
end

M0 = 1;  % equilibrium magnetization
gamma = bmagn/hbar*g; % gyromagnetic ratio

logmsg(1,'Determination of maximum Fourier order (kmax)');
if ~isfield(Options,'kmax') || isempty(Options.kmax)
  % Estimator for maximum relevant Fourier order
  % (1) based on maximum field offset
  maxfieldoffset = max(ModAmp/2-deltaB0,deltaB0+ModAmp/2);
  maxfieldoffset = min(maxfieldoffset,ModAmp);
  maxfreqoffset = (gamma/2/pi)*maxfieldoffset;
  kmax = maxfreqoffset/ModFreq;
  kmax = ceil(kmax*1.4);
  logmsg(1,'  based on max. field offset: %d',kmax);
  
  % (2) based on T2
  %envelope = exp(-2/gam/ModAmp/T2*abs(fourierorder));
  threshold = 1e-6;
  thresholdorder = -log(threshold)/2*gamma*ModAmp*T2;
  thresholdorder = ceil(thresholdorder);
  logmsg(1,'  based on T2: %d',thresholdorder);
  
  % Combine (1) and (2)
  minkmax = 20;
  kmax = min(kmax,thresholdorder);
  kmax = max(kmax,minkmax); % at least 20
  logmsg(1,'  smaller of the two, but at least %d: %d',minkmax,kmax);
  
else
  kmax = Options.kmax; % largest Fourier order, |k|
  logmsg(1,'  user-supplied: %d',kmax);
end


% Solve Bloch equation for steady-state in frequency domain
%-------------------------------------------------------------------------------
% derived by Andrew Ho, December 2013 and January 2014

logmsg(1,'Frequency-domain steady-state solution');
logmsg(1,'  Set up diagonals');
omegam = 2*pi*ModFreq; % modulation angular frequency, rad/s

A = 1i*(-kmax-1:kmax+1).'*omegam;
B = gamma*deltaB0;
C = gamma*B1;
D = gamma*ModAmp/4;
tau1 = 1./(A+1/T1);
tau2 = 1./(A+1/T2);

q = 2:2*kmax+2;
c1 = -D^2*tau2(q-1);
c2 = -B*D*(tau2(q) + tau2(q-1));
c3 = -A(q) - 1/T2 - B^2*tau2(q) - C^2*tau1(q) - D^2*tau2(q-1) - D^2*tau2(q+1);
c4 = -B*D*(tau2(q) + tau2(q+1));
c5 = -D^2*tau2(q+1);
d  = -C*M0*tau1(kmax+2)/T1;

% Assemble pentadiagonal coefficient matrix and RHS vector and solve sparse
% linear system of equations X*Z = Y to get the Fourier coefficients Z = Mky
logmsg(1,'  Assemble sparse pentadiagonal matrix');
nRows = 2*kmax+1;
X = spdiags([c1 c2 c3 c4 c5],[0 1 2 3 4],nRows,nRows+4);
X = X(:,3:end-2);
Y = sparse(nRows,1);
Y(kmax+1) = d;

logmsg(1,'  Solve sparse linear system for y Fourier coefficients');
logmsg(1,'    matrix size: %d x %d',size(X,1),size(X,2));
Z = X\Y;

% Calculate Mkx and Mkz from Mky
logmsg(1,'  Calculate x and z Fourier coefficients');
temp = zeros(2*kmax-1,1);
temp(kmax) = 1;
q = 2:2*kmax;
Mky = Z(q);
if ~onlyAbsorption
  Mkx = tau2(q).*(B*Z(q) + D*Z(q-1) + D*Z(q+1));
  Mkz = tau1(q).*(-C*Z(q) + M0*temp/T1);
end

% Sparse-to-full conversion (since ifft does not support sparse)
Mky = full(Mky);
if ~onlyAbsorption
  Mkx = full(Mkx);
  Mkz = full(Mkz);
end

% Compute time evolution of components
%-------------------------------------------------------------------------------
logmsg(1,'Calculation of time-domain signal.');
tmax = 1/ModFreq; % modulation period

% number of points in time domain
nPoints = Options.nPoints;
if isempty(nPoints)
  nPoints = 2*kmax-1;
end
%nPoints = max(nPoints,2*kmax-1);

t = linspace(0,tmax,nPoints+1).';
t(end) = [];

switch Options.Method
  case 'td' % time-domain explicit construction with complex exponentials
    logmsg(1,'  Method: explicit evolution with complex exponentials');
    if onlyAbsorption
      My = 0;
      for idx = 1:2*kmax-1
        k = idx-kmax; % Fourier order, -(kmax-1):kmax-1
        phase = exp(1i*k*omegam*t);
        My = My + Mky(idx)*phase;
      end
      My = real(My);
    else
      Mx = 0;
      My = 0;
      Mz = 0;
      for idx = 1:2*kmax-1
        k = idx-kmax; % Fourier order, -(kmax-1):kmax-1
        phase = exp(1i*k*omegam*t);
        Mx = Mx + Mkx(idx)*phase;
        My = My + Mky(idx)*phase;
        Mz = Mz + Mkz(idx)*phase;
      end
      Mx = real(Mx);
      My = real(My);
      Mz = real(Mz);
    end
    
  case 'fft' % Computing time domain signals using Fourier transform
    logmsg(1,'  Method: Inverse Fourier transform');
    logmsg(1,'    nPoints: %d, Mky: %d',nPoints,numel(Mky));
    n = numel(Mky);
    if (nPoints<=n)
      % Fourier transform
      My = n*real(ifft(ifftshift(Mky),'symmetric'));
      if ~onlyAbsorption
        Mx = n*real(ifft(ifftshift(Mkx),'symmetric'));
        Mz = n*real(ifft(ifftshift(Mkz),'symmetric'));
      end
      if nPoints<n
        tt = linspace(0,tmax,n+1);
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
      My = FT(Mky);
      if ~onlyAbsorption
        Mx = FT(Mkx);
        Mz = FT(Mkz);
      end
    end
    
  otherwise
    error('Unknown method.');
end

addLastPoint = false;
if addLastPoint
  t(end+1) = 1/ModFreq;
  My(end+1) = My(1);
  if ~onlyAbsorption
    Mx(end+1) = Mx(1);
    Mz(end+1) = Mz(1);
  end
end

t = t/1e-6;    % seconds -> microseconds, for output and plotting

% Graphical rendering
%-----------------------------------------------------------
switch nargout
  case 2
    varargout = {t,My};
  case 4
    varargout = {t,Mx,My,Mz};
  case 7
    varargout = {t,Mx,My,Mz,Mkx,Mky,Mkz};
  case 0
    subplot(2,1,1);
    plot(t,My);
    axis tight;
    xlabel('time (\mus)');
    ylabel('absorption');
    title('Steady-state time-domain signal over one modulation period');
    
    subplot(2,1,2);
    plot(t,Mx);
    axis tight;
    xlabel('time (\mus)');
    ylabel('dispersion');

end

logmsg(1,'=end=blochsteady========%s=================\n',datestr(now));

clear global EasySpinLogLevel
