%fit_levmar  Levenberg-Marquardt nonlinear least squares
%
%   xfit = levmar(funfcn,x0)
%   [xfit,Info] = levmar(funfcn,x0)
%   ... = levmar(funfcn,x0,Opt)
%   ... = levmar(funfcn,x0,Opt,p1,p2,...)
%
%   Find  xm = argmin{F(x)} , where  x = [x_1, ..., x_n]  and
%   F(x) = 0.5 * sum(f_i(x)^2). The functions  f_i(x) (i=1,...,m)
%   must be given by a Matlab function with declaration
%              function  f = funfcn(x, p1,p2,...)
%   p1,p2,... are parameters of the function.  In connection with 
%   nonlinear data fitting they may be arrays with coordinates of 
%   the data points.
%
% Input
%   funfcn  :  Handle to the function.
%   x0   :  Starting guess for xm.
%   Opt  :  options
%         Opt.lambda  starting value of Marquardt parameter
%         Opt.Gradient  termination threshold for gradient
%         Opt.TolStep   termination threshold for step
%         Opt.maxTime   termination threshold for time
%         Opt.delta     step width for Jacobian approximation
%           delta  "relative" step length for difference approximations.
%   p1,p2,..  are passed directly to the function FUN .    
%
% Output
%   xfit :  Computed solution vector.
%   Info :  Performance information, vector with 7 elements:
%           info(1:4) = final values of 
%               [F(x)  ||F'||inf  ||dx||2  mu/max(A(i,i))] ,
%             where  A = Je'* Je .
%           info(5) = no. of iteration steps
%           info(6) = 1 :  Stopped by small gradient
%                     2 :  Stopped by small x-step
%                     3 :  No. of iteration steps exceeded 
%                    -4 :  Dimension mismatch in x, f, B0
%                    -5 :  Overflow during computation
%                    -6 :  Error in approximate Jacobian
%           info(7) = no. of function evaluations

% Method:
% Approximate Gauss-Newton with Levenberg-Marquardt damping and 
% successive updating of Jacobian approximation.

function  [x,info] = fit_levmar(funfcn, x0, Opt, varargin)

if (nargin==0), help(mfilename); return; end
if (nargin<2), error('Need at least 2 arguments!'); end
if (nargin<3),  Opt = []; end

% lambda = starting value of Marquardt parameter
if ~isfield(Opt,'lambda'), Opt.lambda = 1e-3; end
% termation tolerance for gradient (small gradient stops)
if ~isfield(Opt,'Gradient'), Opt.Gradient = 1e-4; end
% termation tolerance for parameter step (small step stops)
if ~isfield(Opt,'TolStep'), Opt.TolStep = 1e-4; end

% delta = relative step for difference approximation
if ~isfield(Opt,'delta'), Opt.delta = 1e-7; end

if ~isfield(Opt,'PrintLevel'), Opt.PrintLevel = 1; end
if ~isfield(Opt,'maxTime'), Opt.maxTime = inf; end


startTime = cputime;

% Check parameters and function call
F = NaN;
norm_g = NaN;
nEvals = 0;

% Check options
delta = Opt.delta;

x0 = x0(:);
n = numel(x0);
if  any(~isreal(x0)) || any(isnan(x0)) || any(isinf(x0)) 
  error('x0 must be real and finite.');
else
  x = x0(:); 
end

stop = 0;

if (~stop)
  [stop F f] = funeval(funfcn,x0,varargin{:});
  nEvals = nEvals + 1;
  if (~stop)
    % Jacobian
    [stop,Je] = JacobianEstimate(funfcn,x,delta,f,varargin{:});
    nEvals = nEvals + n;
    % Check gradient and J'*J
    if (~stop)
      g = Je'*f;
      norm_g = norm(g,inf);
      A = Je'*Je;
      if  isinf(norm_g) || isinf(norm(A(:),inf)), stop = -5; end
    end
  end
end

stop = 0;

if (stop)
  X = x0;
  info.F = F;
  info.norm_g = norm_g;
  info.stop = stop;
  info.nEvals = nEvals;
  return
end

%  Finish initialization
mu = Opt.lambda * max(diag(A)); % initial damping parameter
nu = 2;

norm_h = 0;
j = 0;  % direction of last update

global UserCommand;

iIteration = 0;
hLogLine = findobj('Tag','logLine');

if Opt.PrintLevel
  if isempty(hLogLine)
    fprintf('  Iter   Evals      Function     Gradient          Step\n');
  end
end

while (~stop)
  
  iIteration = iIteration + 1;
  
  if  (norm_g<=Opt.Gradient), stop = 1; break; end
  
  % Levenberg-Marquardt: Compute step and new damping factor
  [h,mu] = ComputeLMStep(A,g,mu);
  norm_h = norm(h);

  if Opt.PrintLevel
    str = sprintf(' %4d:   %5d  %0.5e    %0.5e    %0.5e',iIteration,nEvals,sqrt(F*2),norm_g,norm_h);
    if isempty(hLogLine)
      disp(str)
    else
      set(hLogLine,'String',str);
    end
  end
  
  %if norm_h<=Opt.TolStep, stop = 2; break; end
  if norm_h<=Opt.TolStep*(Opt.TolStep + norm(x)), stop = 2; break; end

  xnew = x + h;
  
  [stop Fnew fnew] = funeval(funfcn,xnew,varargin{:});
  nEvals = nEvals+1;
  if (stop), break; end

  % Update Jacobian estimate Je
  j = mod(j,n) + 1;
  gamma = 0.8;
  if (abs(h(j))<gamma*norm_h)  % recompute with finite differences
    xu = x;
    xu(j) = x(j) + delta;
    [stop Fu fu] = funeval(funfcn,xu,varargin{:});
    nEvals = nEvals+1;
    if (~stop)
      hu = xu - x;
      Je = Je + ((fu-f-Je*hu)/(hu'*hu))*hu';
    end
  end
  Je = Je + ((fnew-f-Je*h)/(h'*h))*h';
  
  % Compute gain ratio
  rho = (F-Fnew)/(0.5*(h'*(mu*h-g)));
  
  % Do step
  if (rho>0)
    x = xnew;
    F = Fnew;
    f = fnew;
  end
  
  % Update damping factor mu
  if (rho>0)
    mu = mu*max(1/3,1-(2*rho-1)^3);
    nu = 2;
  else
    mu = mu*nu;
    nu = 2*nu;
  end
  
  g = Je'*f;
  norm_g = norm(g,inf);
  A = Je'*Je;
  
  if  isinf(norm_g) || isinf(norm(A(:),inf)), stop = -5; break; end
  if (UserCommand==1 || UserCommand==4 || UserCommand==99), stop = 4; break; end
  elapsedTime =  (cputime-startTime)/60;
  if (elapsedTime>Opt.maxTime), stop = 3; break; end

end

if (stop<0)
  Opt.lambda = NaN;
else
  Opt.lambda = mu/max(diag(A));
end

switch (stop)
  case 1, msg = sprintf('Gradient below threshold of %g',Opt.Gradient);
  case 2, msg = sprintf('Parameter step below threshold of %g',Opt.TolStep);
  case 3, msg = sprintf('Time limit of %f minutes reached',Opt.maxTime);
  case 4, msg = sprintf('Stopped by user');
end

if Opt.PrintLevel
  fprintf('Terminated: %s\n',msg);
end

%info = [F norm_g norm_h Opt.lambda iIteration-1 stop nEvals];
info.F = F;
info.norm_g = norm_g;
info.norm_h = norm_h;
info.lambda = Opt.lambda;
info.nIter = iIteration-1;
info.stop = stop;
info.nEvals = nEvals;
info.msg = msg;

return
%======================================================================



%======================================================================
function  [err, J] = JacobianEstimate(funfcn,x,d,f,varargin)
% Approximate Jacobian by finite differences
% Jacobian:
%    dy1/dx1    dy1/dx2   ...
%    dy2/dx1    dy2/dx2   ...
%    ...

nx = numel(x);
J = zeros(numel(f),nx);

for j = 1:nx
  x1 = x;
  x1(j) = x(j)+d;
  fp = feval(funfcn,x1,varargin{:});
  J(:,j) = (fp-f)/d;
end

% Check J
if  ~isreal(J) | any(isnan(J(:))) | any(isinf(J(:)))
  err = -6;
else
  err = 0;
end


%======================================================================
function  [h, mu] = ComputeLMStep(A,g,mu)
% Solve  (Ah + mu*I)h = -g  with possible adjustment of  mu

% Factorize with check of pos. def.
n = size(A,1);  chp = 1;
while  chp
  [R chp] = chol(A + mu*eye(n));
  if  chp == 0  % check for near singularity
    chp = rcond(R) < 1e-15;
  end
  if  chp,  mu = 10*mu; end
end

% Solve  (R'*R)h = -g
h = R \ (R' \ (-g));

%======================================================================
function  [err, F,f] = funeval(funfcn,x,varargin)
%funeval  Check Matlab function which is called by a 
% nonlinear least squares solver.

err = 0;
F = NaN;
n = numel(x);

f = feval(funfcn,x,varargin{:});

sf = size(f);
if  sf(2) ~= 1 | ~isreal(f) | any(isnan(f(:))) | any(isinf(f(:)))
  error('f is not real-valued.');
end

% Objective function
F = (f'*f)/2;
%F = sqrt(mean(f.^2));
if  isinf(F),  err = -5; end
