%fit_simplex    Nelder/Mead downhill simplex minimization algorithm
%
%   xmin = fit_simplex(fcn,x0,FitOpt,...)
%   [xmin,info] = ...
%
%   Tries to minimize fcn(x), starting at x0. FitOpt are options for
%   the minimization algorithm. Any additional parameters are passed
%   to fcn, which can be a string or a function handle.
%
%   Options:
%     Opt.SimplexPars    [rho chi psi sigma], default [1 2 0.5 0.5]
%     Opt.maxTime        maximum time allowed, in minutes
%     Opt.delta              initial extension of simplex

function [x,info] = fit_simplex(errfcn,x0,Opt,varargin)

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
sigma = Opt.SimplexPars(4); % reduction coefficient

if ~isfield(Opt,'Display'), Opt.Display = 'notify'; end
if ~isfield(Opt,'TolStep'), Opt.TolStep = 1e-4; end
if ~isfield(Opt,'TolFun'), Opt.TolFun = 1e-4; end
if ~isfield(Opt,'PrintLevel'), Opt.PrintLevel = 1; end

if ~isfield(Opt,'IterationPrintFunction') || ...
  isempty(Opt.IterationPrintFunction), Opt.IterationPrintFunction = @(str)str; end

%constrain = @(x)sin(x*pi/2); unconstrain = @(x)acos(x)*2/pi;
constrain = @(x)max(min(x,+1),-1); unconstrain = @(x)x;
%constrain = @(x)x; unconstrain = constrain;
n = numel(x0);

% Set up a simplex near the initial guess.
v = repmat(x0(:),1,n+1);
for j = 2:n+1
  v(j-1,j) = v(j-1,j)+delta;
end
%v = v - repmat(mean(v,2),1,n+1); % center

iIteration = 0;
startTime = cputime;

if (Opt.PrintLevel)
  Opt.IterationPrintFunction('initial simplex...');
end

% Evaluate vertices of the simplex
for iVertex = 1:n+1
  x = constrain(v(:,iVertex));
  fv(iVertex) = errfcn(x,varargin{:});
end

% sort so v(1,:) is the best vertex
[fv,idx] = sort(fv); v = v(:,idx);

Procedure = 'initial simplex';
iIteration = iIteration + 1;

if Opt.PrintLevel
  template = ' %4d:  %0.5e  %0.5e  %s';
  str = sprintf(template,iIteration,fv(1),delta,Procedure);
  Opt.IterationPrintFunction(str);
end

% Main algorithm: iterate until
% 1. the maximum distance between the current best vertex and the
%    other vertex in the simplex is less than or equal to TolStep, and
% 2. the corresponding difference in function values is less than or
%    equal to TolFun.

while 1

  % Check whether to stop the loop
  %-----------------------------------------------------------
  elapsedTime = (cputime-startTime)/60;
  if (elapsedTime>Opt.maxTime), stopCode = 1; break; end
  if (UserCommand==1 || UserCommand==4 || UserCommand==99), stopCode = 2; break; end
  if (max(abs(fv(1)-fv(2:n+1))) <= Opt.TolFun) && ...
     (max(max(abs(v(:,2:n+1)-v(:,ones(1,n))))) <= Opt.TolStep)
    stopCode = 3;
    break;
  end

  % Calculate reflection point
  xbar = mean(v(:,1:end-1),2); % average of the n best points
  xr = (1+rho)*xbar - rho*v(:,end);
  xr = constrain(xr);
  fr = errfcn(xr,varargin{:});

  doReduction = 0;
  if fr<fv(:,1) % reflection point is better than best vertex
    % Calculate expansion point
    xe = (1+rho*chi)*xbar - rho*chi*v(:,end);
    xe = constrain(xe);
    fe = errfcn(xe,varargin{:});
    if fe<fr
      Procedure = 'expansion';
      v(:,end) = xe;
      fv(:,end) = fe;
    else
      Procedure = 'reflection';
      v(:,end) = xr;
      fv(:,end) = fr;
    end
  elseif fr<fv(:,n) % reflection point is better than second worst vertex
    Procedure = 'reflection';
    v(:,end) = xr;
    fv(:,end) = fr;
  elseif fr<fv(:,end) % reflection point is better than worst vertex
    % Calculate outside contraction point
    xco = (1+psi*rho)*xbar - psi*rho*v(:,end);
    xco = constrain(xco);
    f = errfcn(xco,varargin{:});
    if f<=fr
      Procedure = 'contraction outside';
      v(:,end) = xco;
      fv(:,end) = f;
    else
      doReduction = 1;
    end
  else % reflection point is worse than (or equal to) worst vertex
    % Calculate inside contraction point
    xci = (1-psi)*xbar + psi*v(:,end);
    xci = constrain(xci);
    f = errfcn(xci,varargin{:});
    if f<fv(:,end)
      Procedure = 'contraction inside';
      v(:,end) = xci;
      fv(:,end) = f;
    else
      doReduction = 1;
    end
  end
 
  if doReduction
    Procedure = 'reduction';
    for iVertex = 2:n+1
      xshr = v(:,1) + sigma*(v(:,iVertex)-v(:,1));
      xshr = constrain(xshr);
      fv(:,iVertex) = errfcn(xshr,varargin{:});
      v(:,iVertex) = xshr;
    end
  end
  
  [fv,idx] = sort(fv);
  v = v(:,idx);
  
  iIteration = iIteration + 1;
  
  if (Opt.PrintLevel)
    thisstep = max(max(abs(v(:,2:n+1)-v(:,ones(1,n)))));
    str = sprintf(template,iIteration,fv(1),thisstep,Procedure);
    Opt.IterationPrintFunction(str);
  end

end

x = v(:,1);
x = unconstrain(x);
info.F = fv(:,1);
info.nIterations = iIteration;
info.elapsedTime = elapsedTime;

switch (stopCode)
  case 1, msg = sprintf('Time limit of %f minutes reached.',Opt.maxTime);
  case 2, msg = sprintf('Stopped by user.');
  case 3, msg = sprintf('Converged (step < %g and function change < %g).',Opt.TolStep,Opt.TolFun);
end

if (Opt.PrintLevel>1)
  fprintf('Terminated: %s\n',msg);
end

return
