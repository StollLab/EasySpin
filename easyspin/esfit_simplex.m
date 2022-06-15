%esfit_simplex    Nelder/Mead simplex minimization algorithm
%
%   xmin = esfit_simplex(fcn,x0,lb,ub,Opt)
%   [xmin,info] = ...
%
%   Tries to find x that minimizes fcn(x), starting at x0. Opt are options
%   for the minimization algorithm. Any additional parameters are passed
%   to fcn, which can be a string or a function handle.
%
%   Input:
%     fcn   ... function to minimize, f(x), where is an array of parameters
%     x0    ... initial parameter values
%     lb    ... lower bounds for parameters
%     ub    ... upper bounds for parameters
%     Opt   ... structure with options
%       .delta         edge length of initial simplex
%       .SimplexPars   [rho, chi, psi, sigma]
%          rho ...     reflection coefficient
%          chi ...     expansion coefficient
%          psi ...     contraction coefficient
%          sigma .     reduction coefficient
%          The default is [1,2,1/2,1/2] for one- and two-dimensional
%          problems, and adaptive for higher dimensions.
%       .maxTime       maximum time allowed, in minutes
%
%   Output:
%     xmin  ... parameter vector with values of best fit
%     info  ... structure with additional information (initial simplex,
%               last simplex, number of iterations, time elapsed)

function [x,info] = esfit_simplex(fcn,x0,lb,ub,Opt)

if nargin==0, help(mfilename); return; end

if nargin<3, Opt = struct; end

% Edge length for initial simplex
if ~isfield(Opt,'delta'), Opt.delta = 0.1; end
delta = Opt.delta;

% Nelder/Mead algorithm parameters
if ~isfield(Opt,'maxTime'), Opt.maxTime = inf; end
if ~isfield(Opt,'SimplexPars')
  Opt.SimplexPars = [1, 2, 0.5, 0.5];
  n = numel(x0);
  if n>2
    % Use adaptive parameters
    % see F. Gao, L. Han, Comput. Optim. Anal. 2010
    % https://doi.org/10.1007/s10589-010-9329-3 
    Opt.SimplexPars = [1, 1+2/n, 0.75-1/(2*n), 1-1/n];
  end
end

rho = Opt.SimplexPars(1); % reflection coefficient
chi = Opt.SimplexPars(2); % expansion coefficient
psi = Opt.SimplexPars(3); % contraction coefficient
sigma = Opt.SimplexPars(4); % reduction coefficient

if ~isfield(Opt,'TolEdgeLength'), Opt.TolEdgeLength = 1e-4; end
if ~isfield(Opt,'TolFun'), Opt.TolFun = 1e-4; end
if ~isfield(Opt,'Verbosity'), Opt.Verbosity = 1; end
if ~isfield(Opt,'IterFcn'), Opt.IterFcn = []; end

if ~isfield(Opt,'IterationPrintFunction') || ...
    isempty(Opt.IterationPrintFunction)
  Opt.IterationPrintFunction = @(str)str;
end

x0 = x0(:);
lb = lb(:);
ub = ub(:);
if numel(lb)~=numel(x0)
  error('Lower-bounds array must have the same number of elements as x0.');
end
if numel(ub)~=numel(x0)
  error('Upper-bounds array must have the same number of elements as x0.');
end
if any(lb>ub)
  error('Lower bounds must not be greater than upper bounds.');
end
nParams = numel(x0);

% Check starting point
if any(~isreal(x0)) || any(isnan(x0)) || any(isinf(x0)) 
  error('x0 must be real and finite.');
end
if any(x0<lb) || any(x0>ub)
  error('Some elements in x0 are out of bounds.');
end

% Transform to (-1,1) interval
transformParams = false;
if transformParams
  transform = @(x) 2*(x-lb)./(ub-lb)-1;
  untransform = @(x) lb + (ub-lb).*(x/2+1/2);
  x0 = transform(x0);
  ub = transform(ub);
  lb = transform(lb);
  fcn = @(x)fcn(untransform(x));
end

constrain = @(x)max(min(x,ub),lb); unconstrain = @(x)x;
%constrain = @(x)sin(x*pi/2); unconstrain = @(x)acos(x)*2/pi;

iIteration = 0;
nEvals = 0;
stopCode = 0;
startTime = cputime;

% Set up a initial simplex near the initial guess.
nVertices = nParams+1;
v = repmat(x0(:),1,nVertices);
v(nParams+1:nParams+1:end) = v(nParams+1:nParams+1:end) + delta.*(ub-lb).';

info.simplex_initial = v;

% Evaluate function at vertices of the simplex
for iVertex = nVertices:-1:1
  x = constrain(v(:,iVertex));
  fv(iVertex) = fcn(x);
  nEvals = nEvals + 1;
end

% Sort so v(1,:) is the best vertex
[fv,idx] = sort(fv);
v = v(:,idx);

Procedure = 'initial simplex';
iIteration = iIteration + 1;

if Opt.Verbosity>0
  logStringFormat = ' iteration %3d: value %0.5e   edge %0.5e   %s';
  str = sprintf(logStringFormat,iIteration,fv(1),delta,Procedure);
  Opt.IterationPrintFunction(str);
end

% Main algorithm: iterate until the following two conditions are both met
% 1. the maximum distance between the current best vertex and the
%    other vertex in the simplex is less than or equal to TolEdgeLength, and
% 2. the corresponding difference in function values is less than or
%    equal to TolFun.

while true
  fv_old = fv;

  % Check whether to stop the iteration loop
  %-----------------------------------------------------------
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>Opt.maxTime
    stopCode = 1;
    break
  end
  largestValDiff = max(abs(fv(1)-fv(2:nParams+1)));
  longestEdgeLength = max(max(abs(v(:,2:nParams+1)-v(:,1))));
  if largestValDiff<=Opt.TolFun && longestEdgeLength<=Opt.TolEdgeLength
    stopCode = 3;
    break
  end
  
  xbar = mean(v(:,1:end-1),2); % average of the n best points
  
  % Calculate reflection point
  xr = xbar + rho*(xbar-v(:,end));
  xr = constrain(xr);
  fr = fcn(xr);
  nEvals = nEvals + 1;

  doReduction = false;
  if fr<fv(:,1) % reflection point is better than best vertex
    % Calculate expansion point
    xe = xbar + rho*chi*(xbar - v(:,end));
    xe = constrain(xe);
    fe = fcn(xe);
    nEvals = nEvals + 1;
    if fe<fr
      Procedure = 'expansion';
      v(:,end) = xe;
      fv(:,end) = fe;
    else
      Procedure = 'reflection';
      v(:,end) = xr;
      fv(:,end) = fr;
    end
  elseif fr<fv(:,nParams) % reflection point is better than second worst vertex
    Procedure = 'reflection';
    v(:,end) = xr;
    fv(:,end) = fr;
  elseif fr<fv(:,end) % reflection point is better than worst vertex
    % Calculate outside contraction point
    %xco = (1+psi*rho)*xbar - psi*rho*v(:,end);
    xco = xbar + psi*rho*(xbar-v(:,end));
    xco = constrain(xco);
    fco = fcn(xco);
    nEvals = nEvals + 1;
    if fco<=fr
      Procedure = 'contraction outside';
      v(:,end) = xco;
      fv(:,end) = fco;
    else
      doReduction = true;
    end
  else % reflection point is worse than (or equal to) worst vertex
    % Calculate inside contraction point
    %xci = (1-psi)*xbar + psi*v(:,end);
    xci = xbar - psi*(xbar-v(:,end));
    xci = constrain(xci);
    fci = fcn(xci);
    nEvals = nEvals + 1;
    if fci<fv(:,end)
      Procedure = 'contraction inside';
      v(:,end) = xci;
      fv(:,end) = fci;
    else
      doReduction = true;
    end
  end
  
  if doReduction
    Procedure = 'reduction';
    for iVertex = 2:nVertices
      xshr = v(:,1) + sigma*(v(:,iVertex)-v(:,1));
      xshr = constrain(xshr);
      fv(:,iVertex) = fcn(xshr);
      nEvals = nEvals + 1;
      v(:,iVertex) = xshr;
    end
  end
  
  % Re-sort vertices according to error values
  [fv,idx] = sort(fv);
  v = v(:,idx);
  
  iIteration = iIteration + 1;
  
  info.bestx = v(:,1);
  info.minF = fv(1);
  info.nEvals = nEvals;
  info.iter = iIteration;
  info.newbest = fv(1)<fv_old(1);
  if ~isempty(Opt.IterFcn)
    UserStop = Opt.IterFcn(info);
  else
    UserStop = false;
  end
  if UserStop
    stopCode = 2;
  end
  
  if Opt.Verbosity>0
    thisstep = max(max(abs(v(:,2:nParams+1)-v(:,ones(1,nParams)))));
    str = sprintf(logStringFormat,iIteration,fv(1),thisstep,Procedure);
    Opt.IterationPrintFunction(str);
  end

  if stopCode~=0, break; end

end

if transformParams
  x = unconstrain(untransform(v(:,1)));
else
  x = unconstrain(v(:,1));
end

info.simplex_final = v;
info.F = fv(:,1);
info.nIterations = iIteration;
info.elapsedTime = elapsedTime;

if Opt.Verbosity>0
  switch stopCode
    case 1, msg = sprintf('Time limit of %f minutes reached.',Opt.maxTime);
    case 2, msg = sprintf('Stopped by user.');
    case 3, msg = sprintf('Converged (edge length < %g and all errors < %g).',Opt.TolEdgeLength,Opt.TolFun);
  end
  fprintf('Terminated: %s\n',msg);
end

info.stop = stopCode;

end
