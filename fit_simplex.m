%fit_simplex    Nelder/Mead downhill simplex minimization algorithm
%
%   xmin = fit_simplex(fcn,x0,Opt,varargin)
%   [xmin,info] = ...
%
%   Tries to minimize fcn(x), starting at x0. Opt are options for
%   the minimization algorithm, any additional parameters are passed
%   to the function to be minimized, fcn, which can be a string or
%   a function handle.
%
%   Options:
%     Opt.SimplexPars    [rho cho psi sigma], default [1 2 0.5 0.5]
%     Opt.maxTime        maximum time allowed, in minutes
%     Opt.delta              initial extension of simplex

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.

function [x,info] = fit_simplex(errfcn,x,Opt,varargin)

if (nargin==0), help(mfilename); return; end

global UserCommand;
if isempty(UserCommand), UserCommand = NaN; end

if (nargin<3), Opt = struct('ununsed',NaN); end

% Parameters for initial simplex
if ~isfield(Opt,'delta'), Opt.delta = 0.1; end
delta = Opt.delta;

% Nelder/Mead algorithm parameters
if ~isfield(Opt,'maxTime'), Opt.maxTime = inf; end
if ~isfield(Opt,'SimplexPars'), Opt.SimplexPars = [1 2 0.5 0.5]; end
rho = Opt.SimplexPars(1); % reflection coefficient
chi = Opt.SimplexPars(2); % expansion coefficient
psi = Opt.SimplexPars(3); % contraction coefficient
sigma = Opt.SimplexPars(4); % shrinkage coefficient

if ~isfield(Opt,'Display'), Opt.Display = 'notify'; end
if ~isfield(Opt,'TolStep'), Opt.TolStep = 1e-4; end
if ~isfield(Opt,'TolFun'), Opt.TolFun = 1e-4; end
if ~isfield(Opt,'PrintLevel'), Opt.PrintLevel = 1; end
tolx = Opt.TolStep;
tolf = Opt.TolFun;

constrain = @(x)sin(x*pi/2); unconstrain = @(x)acos(x)*2/pi;
constrain = @(x)max(min(x,1),-1); unconstrain = @(x)x;
constrain = @(x)x; unconstrain = constrain;
n = numel(x);

onesn = ones(1,n);
two2np1 = 2:n+1;

% Set up a simplex near the initial guess.
v = repmat(x(:),1,n+1);
for j = 2:n+1
  v(j-1,j) = v(j-1,j)+delta;
end
%v = v - repmat(mean(v,2),1,n+1);

iIteration = 0;

hLogLine = findobj('Tag','logLine');
set(hLogLine,'String','initial simplex...');

% Evaluate vertices of the n+1 simplex
for j=1:n+1
  fv(j) = errfcn(v(:,j),varargin{:});
end

% sort so v(1,:) has the lowest function value
[fv,j] = sort(fv); v = v(:,j);

Procedure = 'initial simplex';
iIteration = iIteration + 1;

if Opt.PrintLevel
  template = ' %4d:    %0.5e    %0.5e    %s';
  str = sprintf(template, iIteration, fv(1), delta, Procedure);
  if isempty(hLogLine)
    disp('  Iter      RMS error         Step        Procedure');
    disp(str);
  else
    set(hLogLine,'String',str);
  end
end

% Main algorithm: iterate until
% (a) the maximum coordinate difference between the current best point and the
% other points in the simplex is less than or equal to TolStep. Specifically,
% until max(||v2-v1||,||v3-v1||,...,||v(n+1)-v1||) <= TolStep,
% where ||.|| is the infinity-norm, and v1 holds the
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)

startTime = cputime;

while 1

  % Check whether to stop the loop
  %-----------------------------------------------------------
  elapsedTime =  (cputime-startTime)/60;
  if (elapsedTime>Opt.maxTime), stopCode = 1; break; end
  if (UserCommand==1 || UserCommand==4 || UserCommand==99), stopCode = 2; break; end
  if (max(abs(fv(1)-fv(two2np1))) <= tolf) && ...
     (max(max(abs(v(:,two2np1)-v(:,onesn)))) <= tolx)
    stopCode = 3;
    break;
  end

  % Calculate the reflection point

  % xbar = average of the n (NOT n+1) best points
  xbar = mean(v(:,1:end-1),2);
  xr = (1+rho)*xbar - rho*v(:,end);
  xr = constrain(xr);
  fxr = errfcn(xr,varargin{:});

  if fxr<fv(:,1)
    % Calculate the expansion point
    xe = (1+rho*chi)*xbar - rho*chi*v(:,end);
    xe = constrain(xe);
    fxe = errfcn(xe,varargin{:});
    if fxe<fxr
      Procedure = 'expansion';
      v(:,end) = xe;
      fv(:,end) = fxe;
    else
      Procedure = 'reflection';
      v(:,end) = xr;
      fv(:,end) = fxr;
    end
  else % fxr >= fv(:,1)
    if fxr<fv(:,n)
      Procedure = 'reflection';
      v(:,end) = xr;
      fv(:,end) = fxr;
    else % fxr>=fv(:,n)
      % Perform contraction
      if fxr < fv(:,end)
        xco = (1+psi*rho)*xbar - psi*rho*v(:,end);
        xco = constrain(xco);
        fxco = errfcn(xco,varargin{:});
        if fxco <= fxr
          Procedure = 'contraction outside';
          v(:,end) = xco;
          fv(:,end) = fxco;
        else
          Procedure = 'shrinkage';
          for j=two2np1
            xshr = v(:,1)+sigma*(v(:,j) - v(:,1));
            xshr = constrain(xshr);
            fv(:,j) = errfcn(xshr,varargin{:});
            v(:,j) = xshr;
          end
        end
      else
        xci = (1-psi)*xbar + psi*v(:,end);
        xci = constrain(xci);
        fxci = errfcn(xci,varargin{:});
        if fxci < fv(:,end)
          Procedure = 'contraction inside';
          v(:,end) = xci;
          fv(:,end) = fxci;
        else
          Procedure = 'shrinkage';
          for j=two2np1
            xshr = v(:,1)+sigma*(v(:,j) - v(:,1));
            xshr = constrain(xshr);
            fv(:,j) = errfcn(xshr,varargin{:});
            v(:,j) = xshr;
          end
        end
      end
    end
  end
  [fv,j] = sort(fv);
  v = v(:,j);
  iIteration = iIteration + 1;
  
  if (Opt.PrintLevel)
    thisstep = max(max(abs(v(:,two2np1)-v(:,onesn))));
    str = sprintf(template, iIteration, fv(1), thisstep, Procedure);
    if isempty(hLogLine)
      disp(str);
    else
      set(hLogLine,'String',str);
    end
  end

end

x = v(:,1);
x = unconstrain(x);
info.F = fv(:,1);

switch (stopCode)
  case 1, msg = sprintf('Time limit of %f minutes reached.',Opt.maxTime);
  case 2, msg = sprintf('Stopped by user.');
  case 3, msg = sprintf('Converged (step < %g and function change < %g).',tolx, tolf);
end

if (Opt.PrintLevel>1)
  fprintf('Terminated: %s\n',msg);
end

return
