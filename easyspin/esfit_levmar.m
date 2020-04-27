%esfit_levmar  Levenberg-Marquardt nonlinear least squares fitting
%
%   xfit = esfit_levmar(funfcn,x0)
%   [xfit,Info] = esfit_levmar(funfcn,x0)
%   ... = esfit_levmar(funfcn,x0,Opt)
%   ... = esfit_levmar(funfcn,x0,Opt,p1,p2,...)
%
%   Find  xm = argmin{F(x)} , where  x = [x_1, ..., x_n]  and
%   F(x) = sum(f_i(x)^2)/2. The functions  f_i(x) (i=1,...,m)
%   must be given by a Matlab function with declaration
%              function  f = funfcn(x,p1,p2,...)
%   p1,p2,... are parameters of the function and can be of any type and size.
%
%   The parameter search range is restricted to -1...+1 along each
%   dimension.
%
% Input
%   funfcn Handle to the function.
%   x0     Starting vector in parameter space
%   Opt    Options structure
%            lambda    starting value of Marquardt parameter
%            Gradient  termination threshold for gradient
%            TolStep   termination threshold for step
%            maxTime   termination threshold for time
%            delta     step width for Jacobian approximation
%   p1,p2,... are passed directly to the function funfcn.    
%
% Output
%   xfit    Converged vector in parameter space
%   Info    Performance information, vector with 7 elements:
%           info(1) = final value of F(x)
%           info(2) = final value of ||F'||inf
%           info(3) = final value of ||dx||2
%           info(4) = final value of mu/max(A(i,i)) with A = Je'* Je
%           info(5) = number of iterations
%           info(6) = 1 :  stopped by small gradient
%                     2 :  stopped by small x-step
%                     3 :  max no. of iterations exceeded 
%                    -4 :  dimension mismatch in x, f, B0
%                    -5 :  overflow during computation
%                    -6 :  error in approximate Jacobian
%           info(7) = number of function evaluations

% Method:
% Approximate Gauss-Newton with Levenberg-Marquardt damping and 
% successive updating of Jacobian approximation.

% Search range bounded to -1...+1.

function  [x,info] = esfit_levmar(funfcn, x0, FitOpt, varargin)

if nargin==0, help(mfilename); return; end
if nargin<2, error('Need at least 2 arguments!'); end
if nargin<3,  FitOpt = []; end

% lambda = starting value of Marquardt parameter
if ~isfield(FitOpt,'lambda'), FitOpt.lambda = 1e-3; end
% termation tolerance for gradient (small gradient stops)
if ~isfield(FitOpt,'Gradient'), FitOpt.Gradient = 1e-4; end
% termation tolerance for parameter step (small step stops)
if ~isfield(FitOpt,'TolStep'), FitOpt.TolStep = 1e-4; end

% delta = relative step for difference approximation
if ~isfield(FitOpt,'delta'), FitOpt.delta = 1e-7; end
delta = FitOpt.delta;

if ~isfield(FitOpt,'PrintLevel'), FitOpt.PrintLevel = 1; end
if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'IterationPrintFunction') || ...
    isempty(FitOpt.IterationPrintFunction)
  FitOpt.IterationPrintFunction = @(str)str;
end

startTime = cputime;

% Check parameters and function call
F = NaN;
norm_g = NaN;
nEvals = 0;

nParams = numel(x0);
lb = -ones(nParams,1);
ub = +ones(nParams,1);
if any(lb>ub)
  error('Lower bounds must not be greater than upper bounds.');
end

% Check starting point
x0 = x0(:);
n = numel(x0);
if any(~isreal(x0)) || any(isnan(x0)) || any(isinf(x0)) 
  error('x0 must be real and finite.');
end
if any(x0<lb) || any(x0>ub)
  error('Some elements in x0 are out of bounds.');
end
x = x0(:); 

stopCode = 0;

if ~stopCode
  [stopCode,F,f] = funeval(funfcn,x,varargin{:});
  nEvals = nEvals + 1;
  if ~stopCode
    % Jacobian
    [stopCode,Je] = JacobianEstimate(funfcn,x,f,delta,varargin{:});
    nEvals = nEvals + n;
    % Check gradient and J'*J
    if ~stopCode
      g = Je'*f;
      norm_g = norm(g,inf);
      A = Je'*Je;
      if  isinf(norm_g) || isinf(norm(A(:),inf))
        stopCode = -5;
      end
    end
  end
end

if stopCode
  info.F = F;
  info.norm_g = norm_g;
  info.stop = stopCode;
  info.nEvals = nEvals;
  return
end

%  Finish initialization
mu = FitOpt.lambda * max(diag(A)); % initial damping parameter
nu = 2;

norm_h = 0;
j = 0;  % direction of last update

global UserCommand
if isempty(UserCommand), UserCommand = NaN; end

iIteration = 0;
while ~stopCode
  
  iIteration = iIteration + 1;
  
  if norm_g<=FitOpt.Gradient
    stopCode = 1;
    break;
  end
  
  % Levenberg-Marquardt: Compute step and new damping factor
  [h,mu] = ComputeLMStep(A,g,mu);
  norm_h = norm(h);

  if FitOpt.PrintLevel
    str = sprintf(' iteration %d:\n    value %0.5f\n    gradient %0.5f\n    step %0.5f',iIteration,sqrt(F*2),norm_g,norm_h);
    FitOpt.IterationPrintFunction(str);
  end
  
  if norm_h<=FitOpt.TolStep*(FitOpt.TolStep + norm(x))
    stopCode = 2;
    break;
  end
  
  xnew = x + h;
  xnew = min(max(xnew,lb),ub); % apply bounds
  
  [stopCode,Fnew,fnew] = funeval(funfcn,xnew,varargin{:});
  nEvals = nEvals+1;
  if stopCode, break; end

  % Update Jacobian estimate Je
  j = mod(j,n) + 1;
  gamma = 0.8;
  if abs(h(j))<gamma*norm_h  % recompute with finite differences
    xu = x;
    xu(j) = x(j) + delta;
    [stopCode,~,fu] = funeval(funfcn,xu,varargin{:});
    nEvals = nEvals+1;
    if ~stopCode
      hu = xu - x;
      Je = Je + ((fu-f-Je*hu)/(hu'*hu))*hu';
    end
  end
  Je = Je + ((fnew-f-Je*h)/(h'*h))*h';
  
  % Compute gain ratio
  rho = (F-Fnew)/(0.5*(h'*(mu*h-g)));
  
  % Do step
  if rho>0
    x = xnew;
    F = Fnew;
    f = fnew;
  end
  
  % Update damping factor mu
  if rho>0
    mu = mu*max(1/3,1-(2*rho-1)^3);
    nu = 2;
  else
    mu = mu*nu;
    nu = 2*nu;
  end
  
  g = Je'*f;
  norm_g = norm(g,inf);
  A = Je'*Je;
  
  if isinf(norm_g) || isinf(norm(A(:),inf)), stopCode = -5; break; end
  if UserCommand==1 || UserCommand==4 || UserCommand==99, stopCode = 4; break; end
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>FitOpt.maxTime, stopCode = 3; break; end

end

if stopCode<0
  FitOpt.lambda = NaN;
else
  FitOpt.lambda = mu/max(diag(A));
end

switch stopCode
  case 1, msg = sprintf('Gradient below threshold of %g',FitOpt.Gradient);
  case 2, msg = sprintf('Parameter step below threshold of %g',FitOpt.TolStep);
  case 3, msg = sprintf('Time limit of %f minutes reached',FitOpt.maxTime);
  case 4, msg = sprintf('Stopped by user');
end

if FitOpt.PrintLevel>1
  fprintf('Terminated: %s\n',msg);
end

info.F = F;
info.norm_g = norm_g;
info.norm_h = norm_h;
info.lambda = FitOpt.lambda;
info.nIter = iIteration-1;
info.stop = stopCode;
info.nEvals = nEvals;
info.msg = msg;

return
%======================================================================



%======================================================================
function  [err, J] = JacobianEstimate(funfcn,x0,f0,delta,varargin)
% Compute approximate Jacobian using finite differences
% Jacobian:
%    dy1/dx1    dy1/dx2   ...
%    dy2/dx1    dy2/dx2   ...
%    ...

nVariables = numel(x0);
J = zeros(numel(f0),nVariables);

for ix = 1:nVariables
  x1 = x0;
  x1(ix) = x0(ix) + delta;
  f1 = funfcn(x1,varargin{:});
  f1 = f1(:);
  J(:,ix) = (f1-f0)/delta;
end

% Check J
if  ~any(isreal(J(:))) || any(isnan(J(:))) || any(isinf(J(:)))
  err = -6;
else
  err = 0;
end


%======================================================================
function [h,mu] = ComputeLMStep(A,g,mu)
% Solve (A+mu*1)*h = -g, scaling mu if needed; using Cholesky factorization

notPosDef = true;
while notPosDef
  [R,notPosDef] = chol(A + mu*eye(size(A)));
  if ~notPosDef
    % check whether close to singular
    notPosDef = rcond(R)<1e-15;
  end
  if notPosDef
    mu = 10*mu;
  end
end

% Solve  (R'*R)*h = -g
h = R\(R'\(-g));

%======================================================================
function  [errCode,F,f] = funeval(funfcn,x,varargin)
%funeval  Check Matlab function which is called by a nonlinear least squares solver.

errCode = 0;

f = funfcn(x,varargin{:});
f = f(:);

if any(isnan(f))
  error('f contains at least one NaN value.')
end

if any(isinf(f))
  error('f contains at least one inf value.')
end

if  any(~isreal(f))
  error('f is not real-valued.');
end

% Objective function
F = f'*f;

if isinf(F), errCode = -5; end
