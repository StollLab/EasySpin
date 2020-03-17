% propint  Compute propagation operator 
%
%   U = propint(H0,H1,t,freq)
%   U = propint(H0,H1,t,freq,phase)
%   U = propint(H0,H1,t,freq,phase,n)
%   [U,SaveU] = propint(...)
%
%   U = propint(SaveU,t,freq)
%   U = propint(SaveU,t,freq,phase)
%
%   Computes a propagator by integrating
%
%      expm(-i*2*pi*(H0+H1*cos(2*pi*freq*t+phase))*t)
%
%   over the interval t(1) to t(2) using steps of
%   dt = 1/freq/n. Defaults: phase = 0, n = 256.
%   t = t2 is equivalent to t = [0 t2]. H0 and
%   H1 should be in frequency unit, t in complementary
%   time units, e.g. GHz and ns or MHz and µs.
%
%   SaveU is a cell array containing intermediate
%   results that can be reused by supplying SaveU
%   instead of H0 and H1 in the next call. The
%   speed-up is huge.
%
%   Example: An resonant pi/2 pulse for a spin-1 would be
%
%      Sz = sop(1,'z'); Sy = sop(1,'y'); % S=1 system
%      mwFreq = 10e3; tp = 0.010; % 10 GHz, 10 ns
%      U = propint(mwFreq*Sz,1/2/tp*Sy,tp,mwFreq)
%
%   Example: A resonance pi pulse on an two spin system
%
%      spins = [1/2 1/2]
%      Sz = sop(spins,'ze'); Sy = sop(spins,'ye');
%      mwFreq = 10e3; tp = 0.010; % 10 GHz, 10 ns
%      U = propint(mwFreq*Sz,1/tp*Sy,tp,mwFreq);
%      Density = -Sz; U*Density*U'

function varargout = propint(varargin)

if (nargin==0), help(mfilename); return; end

% log messages for developers...
Display = 0;

%Method = 0; % simple method, used 
%Method = 1; % Runge-Kutta 4, Oct 2007, Alexey Potapov
Method = 0;

% argument parsing
%----------------------------------------------------------------
switch nargout
case 0, SaveU = 0;
case 1, SaveU = 0;
case 2, SaveU = 1;
otherwise
  error('Wrong number of output arguments');
end

UseUstore = iscell(varargin{1});
if UseUstore
  [Ustore,tlim,Frequency] = deal(varargin{1:3});
  nIntervals = numel(Ustore)-1;
  if nargin>3, Phase = varargin{4}; else Phase = 0; end
else
  [H0,H1,tlim,Frequency] = deal(varargin{1:4});
  if nargin>4, Phase = varargin{5}; else Phase = 0; end
  if nargin>5, nIntervals = varargin{6}; else nIntervals = 256; end
end

if numel(tlim)==1, tlim = [0 tlim]; end

if (numel(Frequency)~=1) | ~isreal(Frequency) | (Frequency<=0)
  error('Frequency must be real and positive!');
end

if (numel(tlim)~=2) | ~isreal(tlim)
  error('tlim is not a real 1x2 array!');
end

% Short-cut if no time dependence
%----------------------------------------------------------------
if ~UseUstore & all(H1(:)==0)
  U = expm(-i*2*pi*diff(tlim)*H0);
  return
end

% Pre-calculations of times and periods...
%----------------------------------------------------------------

tPeriod = 1/Frequency; % length of one period
dt = tPeriod/nIntervals; % length of one integration interval
t = (0.5:nIntervals)*dt; % time points at the center of the intervals
t = (0.5:nIntervals)*dt; % time points at the center of the intervals

nSteps = diff(tlim)/dt; % #integration step theoretically neede
if abs(fix(nSteps)/nSteps-1)>1e-2
  % Pulse length gets rounded to next integer nStep. If this is
  % a serious change, issue warning!
  warning('Pulse too short. Increase number of integration points!');
end

% Indices to intervals where start and end of pulse fall.
j = fix(nIntervals*rem(tlim/tPeriod,1)+1);
nPeriods = tlim/tPeriod; % number of periods the pulse lasts.

% Do some log display.
if Display
  fprintf('frequency %0f GHz, one period is %0g ns\n',...
	  Frequency*1e-3,tPeriod*1e3);
  fprintf('pulse length %0g ns\n',diff(tlim)*1e3);
  fprintf('  is %0f integration steps of %0g ns\n',nSteps,dt*1e3);
  fprintf('start t %0g ns, %0g periods over, j %d\n',...
	  tlim(1)*1e3,nPeriods(1),j(1)-1);
  fprintf('end   t %0g ns, %0g periods over, j %d\n',...
	  tlim(2)*1e3,nPeriods(2),j(2)-1);
  fprintf('periods %0g\n',diff(nPeriods));
  fprintf('  (%0g)(%d)(-%0g)',...
	  rem(nPeriods(2),1),diff(fix(nPeriods)),rem(nPeriods(1),1));
end


% Integrate propagator over one period
%----------------------------------------------
if ~UseUstore

  % Start propagator U(0,0) is identity matrix
  UZero = eye(size(H0));
  UPeriod = UZero;

  % Pre-compute operators incl. constants, cos function
  cH0 = -i*2*pi*dt*H0;
  cH1 = -i*2*pi*dt*H1;
  switch Method
    case 0
      ct = cos(2*pi*Frequency*t+Phase);
    case 1
      ct_1 = cos(2*pi*Frequency*(t-dt/2)+Phase); % k1
      ct  =  cos(2*pi*Frequency*t+Phase);        % k2==k3
      ct_2 = cos(2*pi*Frequency*(t+dt/2)+Phase); %  k4
  end

  % Allocate space for intermediate propagator results
  if SaveU
    % n intervals give n+1 propagators (defined at the borders)
    Ustore = cell(1,nIntervals+1);
    % The first propagator is the identity matrix.
    Ustore{1} = UPeriod;
  end

  % Loop over all intervals.
  for k = 1:nIntervals
    % Save propagators for start and end partial period
    if k==j(1), U1 = UPeriod; end
    if k==j(2), U2 = UPeriod; end

    % Integrate over next interval
    switch Method
      case 0
        UPeriod = expm(cH0 + cH1*ct(k)) * UPeriod;
      case 1
        k1 = ( cH0 + cH1*ct_1(k) )*UPeriod;
        k2 = ( cH0 + cH1*ct(k)   )*(UPeriod+k1/2);
        k3 = ( cH0 + cH1*ct(k)   )*(UPeriod+k1/2);
        k4 = ( cH0 + cH1*ct_2(k) )*(UPeriod+k3);
        UPeriod =  UPeriod +1/6*(k1+2*k2+2*k3+k4);
    end
    
    % Store if wanted.
    if SaveU, Ustore{k+1} = UPeriod; end
  end
else
  % Retrieve propagators needed for computation.
  UPeriod = Ustore{end};
  U1 = Ustore{j(1)+1};
  U2 = Ustore{j(2)+1};
end

% Computation of total propagator U
%----------------------------------------------------------------
% (a) U = U1'
%     back-propagate to period border prior to tlim(1)
% (b) U = UPeriod^diff(fix(nPeriods)) * U
%     Propagate by diff(fix(nPeriods)) periods to
%     period border prior to tlim(2)
% (c) U = U2 * U
%     Propagate to tlim(2)

U  = U2 * UPeriod^diff(fix(nPeriods)) * U1';


% Arrange output
%----------------------------------------------------------------
switch nargout
case 1, varargout = {U};
case 2, varargout = {U,Ustore};
end

return
