% propint  Compute propagation matrix
%
%   U = propint(H0,H1,tlim,freq)
%   U = propint(H0,H1,tlim,freq,phase)
%   U = propint(H0,H1,tlim,freq,phase,Opt)
%   [U,Ustore] = propint(...)
%
%   U = propint(Ustore,tlim,freq)
%   U = propint(Ustore,tlim,freq,phase)
%
%   Computes a propagator matrix U by integrating
%
%      dU/dt = -2i*pi*(H0+H1*cos(2*pi*freq*t+phase)*U
%
%   over the interval t = tlim(1) to tlim(2).
%
%   H0 and H1 are Hermitian matrices and should be in MHz.
%
%   tlim is a 2-element array [t1 t2] indicating start and end time of the 
%   integration, in µs. Alternatively, tlim = t2 is equivalent to tlim = [0 t2].
%
%   freq and phase are the frequency and the phase of the oscillating part of
%   the Hamiltonian. freq is in MHz, and phase in radians.
%
%   Opt specifies settings for the integration algorithm:
%   .Method   integration method ("simple" or "RK4") (default "simple")
%             "simple" uses the trapezoidal rule
%             "RK4" uses the 4th-order Runge-Kutta method
%   .n        number of intervals over one period (default 360)
%
%   Ustore is a cell array containing intermediate results that can be
%   reused by supplying Ustore instead of H0 and H1 in the next call. The
%   speed-up is massive.
%
%   Example: A resonant pi/2 pulse for a spin-1/2 would be
%
%      [Sz,Sy] = sop(1/2,'z','y');  % spin operator matrices
%      nu0 = 10e3;           % resonance frequency, MHz
%      tp = 0.010;           % pulse duration, µs
%      theta = pi/2;         % pulse flip angle, rad
%      nu1 = theta/pi/tp;    % pulse amplitude, MHz
%      nu_mw = 10e3;         % microwave frequency, MHz
%      H0 = nu0*Sz;          % static Hamiltonian, MHz
%      H1 = nu1*Sy;          % microwave Hamiltonian MHz
%      U = propint(H0,H1,tp,nu_mw)

% Integration method:
%  "simple"  trapezoidal rule (default)
%  "RK4"     Runge-Kutta 4th order (added by Alexey Potapov, Oct 2007)

function varargout = propint(varargin)

if nargin==0
  help(mfilename);
  return
end

% log messages for developers
diagnostics = false;

% Argument parsing
%-------------------------------------------------------------------------------
switch nargout
  case 0, storeU = false;
  case 1, storeU = false;
  case 2, storeU = true;
otherwise
  error('Wrong number of output arguments');
end

useUstore = iscell(varargin{1});
if useUstore
  [Ustore,tlim,frequency] = deal(varargin{1:3});
  if nargin>3, phase = varargin{4}; else, phase = 0; end
  Opt = struct;
  Opt.n = numel(Ustore)-1;
else
  if nargin<4
    error('At least 4 inputs are required (H0, H1, tlim, freq).');
  end
  [H0,H1,tlim,frequency] = deal(varargin{1:4});
  phase = 0;
  Opt = struct;
  if nargin>=5, phase = varargin{5}; end
  if isstruct(phase)
    error('5th input (phase) must be a scalar, in radians. You provided a structure.');
  end
  if nargin==6, Opt = varargin{6}; end
  if nargin>6
    error('No more than 6 input arguments are possible.');
  end
end

% Check time range input
if isscalar(tlim)
  tlim = [0 tlim];
end
if numel(tlim)~=2 || ~isreal(tlim) || diff(tlim)<=0
  error('t is must be a real-valued scalar or a real-valued 2-element array with t(1)<t(2).');
end

% Check frequency input
if numel(frequency)~=1 || ~isreal(frequency) || frequency<=0
  error('Frequency must be real and positive.');
end

% Number of points per period
if ~isfield(Opt,'n') || isempty(Opt.n)
  Opt.n = 360;
end
nIntervals = Opt.n;

if ~isfield(Opt,'Method')
  Opt.Method = "simple";
end

% Short-cut if no time dependence
%-------------------------------------------------------------------------------
if ~useUstore && all(H1(:)==0)
  U = expm(-2i*pi*diff(tlim)*H0);
  switch nargout
    case {0,1}
      varargout = {U};
    case 2
      varargout = {U,Ustore};
  end
  return
end

% Pre-calculations of times and periods
%-------------------------------------------------------------------------------
tPeriod = 1/frequency;  % length of one period
dt = tPeriod/nIntervals;  % length of one integration interval
t = (0.5:nIntervals)*dt;  % time points at the center of the intervals

nSteps = diff(tlim)/dt; % numer of integration step theoretically needed
if abs(fix(nSteps)/nSteps-1)>1e-2
  % Pulse length gets rounded to next integer nStep. If this is
  % a serious change, issue warning!
  warning('Pulse too short. Increase number of integration points!');
end

% Indices to intervals where start and end of pulse fall.
j = fix(nIntervals*rem(tlim/tPeriod,1)+1);
nPeriods = tlim/tPeriod; % number of periods the pulse lasts.

% Print diagnostic information
if diagnostics
  fprintf('frequency %0f GHz, one period is %0g ns\n',...
	  frequency*1e-3,tPeriod*1e3);
  fprintf('pulse length %0g ns\n',diff(tlim)*1e3);
  fprintf('  = %f integration steps of %g ns\n',nSteps,dt*1e3);
  fprintf('start t %0g ns, %0g periods over, j %d\n',...
	  tlim(1)*1e3,nPeriods(1),j(1)-1);
  fprintf('end   t %0g ns, %0g periods over, j %d\n',...
	  tlim(2)*1e3,nPeriods(2),j(2)-1);
  fprintf('periods %0g\n',diff(nPeriods));
  fprintf('  (%0g)(%d)(-%0g)\n',...
	  rem(nPeriods(2),1),diff(fix(nPeriods)),rem(nPeriods(1),1));
end


% Integrate propagator over one period
%-------------------------------------------------------------------------------
if ~useUstore

  % Start propagator U(0,0) is identity matrix
  Uzero = eye(size(H0));
  Uperiod = Uzero;

  % Pre-compute matrices including constants and cos function
  cH0 = -2i*pi*dt*H0;
  cH1 = -2i*pi*dt*H1;
  switch Opt.Method
    case "simple"
      ct = cos(2*pi*frequency*t+phase);
    case "RK4"
      ct_1 = cos(2*pi*frequency*(t-dt/2)+phase); % k1
      ct  =  cos(2*pi*frequency*t+phase);        % k2==k3
      ct_2 = cos(2*pi*frequency*(t+dt/2)+phase); % k4
  end

  % Allocate space for intermediate propagator results
  if storeU
    % n intervals give n+1 propagators (defined at the borders)
    Ustore = cell(1,nIntervals+1);
    % The first propagator is the identity matrix.
    Ustore{1} = Uperiod;
  end

  % Loop over all intervals
  for k = 1:nIntervals
    % Save propagators for start and end partial period
    if k==j(1), U1 = Uperiod; end
    if k==j(2), U2 = Uperiod; end

    % Integrate over next interval
    switch Opt.Method
      case "simple"
        Uperiod = expm(cH0 + cH1*ct(k)) * Uperiod;
      case "RK4"
        k1 = (cH0 + cH1*ct_1(k))*Uperiod;
        k2 = (cH0 + cH1*ct(k)  )*(Uperiod+k1/2);
        k3 = (cH0 + cH1*ct(k)  )*(Uperiod+k1/2);
        k4 = (cH0 + cH1*ct_2(k))*(Uperiod+k3);
        Uperiod =  Uperiod + (k1+2*k2+2*k3+k4)/6;
    end
    
    % Store if wanted
    if storeU, Ustore{k+1} = Uperiod; end
  end

else
  
  % Retrieve propagators needed for computation
  Uperiod = Ustore{end};
  U1 = Ustore{j(1)+1};
  U2 = Ustore{j(2)+1};
  
end

% Computation of total propagator U
%-------------------------------------------------------------------------------
% (a) U = U1'
%     back-propagate to period border prior to tlim(1)
% (b) U = UPeriod^diff(fix(nPeriods)) * U
%     Propagate by diff(fix(nPeriods)) periods to
%     period border prior to tlim(2)
% (c) U = U2 * U
%     Propagate to tlim(2)

U  = U2 * Uperiod^diff(fix(nPeriods)) * U1';


% Arrange output
%----------------------------------------------------------------
switch nargout
  case {0,1}
    varargout = {U};
  case 2
    varargout = {U,Ustore};
end

end
