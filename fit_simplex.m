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

function [x,info] = fit_simplex(funfcn,x,Opt,varargin)

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

use_feval = sscanf(version,'%s',1)<7;

% Convert to function handle as needed.
funfcn = fcnchk(funfcn,length(varargin));
n = numel(x);

onesn = ones(1,n);
two2np1 = 2:n+1;
one2n = 1:n;

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = zeros(n,n+1);
fv = zeros(1,n+1);
v(:,1) = xin;
x(:) = xin;    % Change x to the form expected by funfcn
fv(:,1) = feval(funfcn,x,varargin{:});
iIteration = 0;

hLogLine = findobj('Tag','logLine');


% Continue setting up the initial simplex.
for j = 1:n
  y = xin;
  y(j) = y(j) + delta;
  v(:,j+1) = y;
  x(:) = y;
  fv(1,j+1) = feval(funfcn,x,varargin{:});
end

% sort so v(1,:) has the lowest function value
[fv,j] = sort(fv);
v = v(:,j);

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
% until max(||v2-v1||,||v2-v1||,...,||v(n+1)-v1||) <= TolStep,
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
  xbar = sum(v(:,one2n), 2)/n;
  xr = (1 + rho)*xbar - rho*v(:,end);
  x(:) = xr; fxr = feval(funfcn,x,varargin{:});

  if fxr < fv(:,1)
    % Calculate the expansion point
    xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
    x(:) = xe;
    fxe = feval(funfcn,x,varargin{:});
    if fxe < fxr
      v(:,end) = xe;
      fv(:,end) = fxe;
      Procedure = 'expansion';
    else
      v(:,end) = xr;
      fv(:,end) = fxr;
      Procedure = 'reflection';
    end
  else % fxr >= fv(:,1)
    if fxr < fv(:,n)
      v(:,end) = xr;
      fv(:,end) = fxr;
      Procedure = 'reflection';
    else % fxr >= fv(:,n)
      % Perform contraction
      if fxr < fv(:,end)
        % Perform an outside contraction
        xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
        x(:) = xc;
        fxc = feval(funfcn,x,varargin{:});

        if fxc <= fxr
          v(:,end) = xc;
          fv(:,end) = fxc;
          Procedure = 'contraction outside';
        else
          % perform a shrink
          Procedure = 'shrinkage';
        end
      else
        % Perform an inside contraction
        xcc = (1-psi)*xbar + psi*v(:,end);
        x(:) = xcc; 
        fxcc = feval(funfcn,x,varargin{:});

        if fxcc < fv(:,end)
          v(:,end) = xcc;
          fv(:,end) = fxcc;
          Procedure = 'contraction inside';
        else
          % perform a shrink
          Procedure = 'shrinkage';
        end
      end
      if strcmp(Procedure,'shrinkage')
        for j=two2np1
          v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
          x(:) = v(:,j);
          fv(:,j) = feval(funfcn,x,varargin{:});
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

x(:) = v(:,1);
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
