% exponfit  Fits exponential(s) to 1D and 2D data
%
%   [k,c] = exponfit(x,y)
%   ...   = exponfit(x,y,nExps)
%   ...   = exponfit(x,y,nExps,'noconst')
%   [k,c,yFit] = exponfit(...)
%   k = exponfit(...)
%
%   Computes least square fits for single or double exponential
%   decays or recoveries. For matrices it works along columns and
%   returns the result in corresponding matrix form. The fitting model
%   is of the form
%
%      c(1) + c(2)*exp(-k(1)*x) + c(3)*exp(-k(2)*x)
%
%   Input:
%     x          abscissa values, a vector
%     y          ordinate values, a vector or a matrix
%     nExps      number of exponents to include (1 or 2, default is 1)
%     'noconst'  if given, c(1) is set to 0 and not included in the fit
%
%   Output:
%     k       exponents
%     c       linear coefficients
%     yFit    the fitting model function

function varargout = exponfit(x,yData,nExponents,noconst)

if (nargin==0), help(mfilename); return; end

if (nargin<2) || (nargin>4), error('Wrong number of input arguments!'); end
if (nargout<1), error('Not enough output arguments.'); end
if (nargout>3), error('Too many output arguments.'); end

if (nargin<3), nExponents = 1; end
if (nargin<4), noconst = '-'; end

switch nExponents
  case {1,2}
  otherwise
    error('Number of exponents must be 1 or 2.');
end

ComputeFitFunc = (nargout>2);

IncludeOffset=  ~strcmp(noconst,'noconst');

Options = [IncludeOffset 1e-11 1e-11 1e4 1e-3];

nRows = size(yData,1);
Flipped = (nRows==1);
if Flipped, yData = yData.'; end
[nRows,nCols] = size(yData);
x = x(:);

% Rescale horizontal axis
scale = max(abs(x));
x = x/scale;

% Pre-allocation.
k = zeros(nExponents,nCols);
c = zeros(nExponents+1,nCols);

if ComputeFitFunc
  yFitted = zeros(nRows,nCols);
else
  yFitted = [];
end

for iCol = 1:nCols
  
  % Get column data and invert if it is not a decay.
  ySlice = yData(:,iCol);
  IsRecovery = ySlice(1)<ySlice(end);
  if IsRecovery, ySlice = -ySlice; end
  
  % Get an initial guess for the decay constant(s).
  kStart = -GuessSingleDecay(x,ySlice);
  if (nExponents==2), kStart = [.8;1.3]*kStart; end
  
  % Do the separable nonlinear least squares fitting.
  [k(:,iCol),c_,info] = mexpfit2(x,ySlice,kStart,Options);
  if ~IncludeOffset
    c(1,iCol) = 0;
    c(2:end,iCol) = c_;
  else
    c(:,iCol) = c_;
  end
  
  % Undo the inversion.
  if IsRecovery
    c(:,iCol) = -c(:,iCol);
  end
  
  
  % Compute fitted function.
  if ComputeFitFunc
    switch nExponents
    case 1
      yFitted(:,iCol) = c(1,iCol) + c(2,iCol)*exp(k(1,iCol)*x);
    case 2
      yFitted(:,iCol) = c(1,iCol) + c(2,iCol)*exp(k(1,iCol)*x) + ...
        c(3,iCol)*exp(k(2,iCol)*x);
    end
  end
  
end

% Rescale
k = k/scale;

if Flipped
  yFitted = yFitted.';
else
  k = k.';
  c = c.';
end

varargout = {-k,c,yFitted};
varargout = varargout(1:nargout);

return

function kGuessed = GuessSingleDecay(t,y)
% estimate coefficients for y = A*exp(-k*x)+Offset
C = y(end);
A = y(1)-C;
dt = t(2)-t(1);
integral_ = sum(y-C)*dt;
kGuessed = A/integral_;
if (kGuessed<0)
  kGuessed = 1/(length(t)/2*dt);
end
return

function  varargout = mexpfit2(tData,yData,xGuess,Options)
% MEXPFIT   Least squares fit to sum of decaying exponentials
%
%   [k,c,info{,perf}] = mexpfit(tData,yData,k0,Options)
%
%   Input
%     tData     abscissa values of data
%     yData     ordinate values of data
%     k0    initial guess for the decay rates
%     Options  options, 4- or 5-element vector
%           Options(1)  0 omit constant term
%           Options(2)  gradient limit
%           Options(3)  step limit
%           Options(4)  no of functions evaluations
%           Options(5)  starting value for Marquardt (1e-3 default)
%
%   Output
%     k     decay rates for least squares fit
%     c     linear coefficients of least squares fit
%     info  information vector
%           info(1:3) = final values of  F(x)  |F'|inf  |dx|2  
%           info(4) = no. of evaluations of (F,Jacobian)
%           info(5) = 1 :  Stopped by small gradient
%                     2 :  Stopped by small x-step
%                     3 :  Stopped by  MaxFunEvals
%                     4 :  Stopped by extreme step
%                     5 :  Stopped by stalling.
%     perf  (optional). If present, then array, holding 
%           perf(1:2,:) = values of  F(x) and || F'(x) ||inf
%           perf(3,:) = mu-values.

% Hans Bruun Nielsen,  IMM, DTU.  00.02.14
% modified by Stefan Stoll, ETH Zurich, 04-apr-01

if (nargout<2), error('Not enough output arguments.'); end
if (nargout>4), error('Too many output arguments.'); end
if (nargin<4) || (nargin>4), error('Wrong number of input arguments!'); end

if numel(tData)~=numel(yData)
  error('tData and yData must have the same number of elements!');
end

%  Initialise
tData = tData(:); % assure column format
yData = yData(:); % assure column format
x = -sort(-xGuess(:)); % assure column format, sort in descending order
n = length(x);

if  any(x >= 0), error('k must be strictly negative'), end
if  length(Options) < 5,  Options(5) = 1e-3; end
if  Options(5) <= 0,  Options(5) = 1e-3; end

IncludeOffset = Options(1);
MinGradNorm = Options(2);
MaxFunEvals = Options(4);

%  Initial values
[Differences,Jacobian,c] = func(x,tData,yData,IncludeOffset);
A = Jacobian'*Jacobian;
Gradient = Jacobian'*Differences;
F = (Differences'*Differences)/2;
GradientNorm = norm(Gradient,inf);
mu = Options(5) * max(diag(A));
Trace = nargout>3;
if  Trace
  X = x(:,ones(1,MaxFunEvals+1));
  Performance = [F; GradientNorm; mu]*ones(1,MaxFunEvals+1);
else
  Performance = [];
end 
iIteration = 0;
nFunEvals = 1;
nu = 2;
StopCriterion = 0;
mok = 0;
NormStep = NaN;

% Iterations: Loop as long as no stop criterion is fulfilled.
while ~StopCriterion
  
  % First check for various stop criteria and compute step
  if  GradientNorm <= MinGradNorm
    StopCriterion = 1;
  elseif  nFunEvals == MaxFunEvals
    StopCriterion = 3;
  else
    Step = (A + mu*eye(n))\(-Gradient);
    NormStep = norm(Step);
    nx = Options(3) + norm(x);
    if  NormStep <= Options(3)*nx
      StopCriterion = 2;
    elseif NormStep >= nx/eps
      StopCriterion = 4;
    end    % Almost singular ?
  end
  
  % If none of the stop criteria is fulfilled, do the next iteration
  if  ~StopCriterion
    iIteration = iIteration + 1;
    %  Check Step
    i = find(x + Step >= 0);
    if  ~isempty(i)    % Reduce Step
      Step = min(-.99*x(i)./Step(i))*Step;
      mu = 10*mu;
    end
    xnew = x + Step;
    Step = xnew - x;
    dL = (Step'*(mu*Step - Gradient))/2; 
    [NewDifferences,NewJacobian] = func(xnew,tData,yData,IncludeOffset);
    nFunEvals = nFunEvals + 1;   
    newF = (NewDifferences'*NewDifferences)/2;
    dF = F - newF;
    if  (dL > 0) & (dF > 0)
      mok = mok + 1;
      mu = mu * max(1/3, 1 - (2*dF/dL - 1)^3);   nu = 2;
      x = xnew;   F = newF;
      Jacobian = NewJacobian;
      Differences = NewDifferences; 
      A = Jacobian'*Jacobian;
      Gradient = Jacobian'*Differences;
      GradientNorm = norm(Gradient,inf);
    else  % Marquardt fail
      mu = mu*nu;
      nu = 2*nu;
      if  (mok>n) & (nu>8)
        StopCriterion = 5;
      end
    end
    if Trace
      X(:,iIteration) = x(:);
      Performance(:,iIteration) = [F; GradientNorm; mu];
    end
  end   
end

%  Set return values
x = -sort(-x);
if  Trace
  X = [X(:,1:iIteration-1) x];
  Performance = Performance(:,1:iIteration);
else
  X = x;
end
[r,Jacobian,c] = func(x, tData, yData, IncludeOffset);
F = .5*dot(r,r);   
if IncludeOffset
  c = c([n+1 1:n]);
end

info = [F,GradientNorm,NormStep,nFunEvals,StopCriterion];
varargout = {X,c,info,Performance};
varargout = varargout(1:nargout);

return


% ==========  auxiliary function  =================================
function  [Residuals,Jacobian,c] = func(x,tData,yData,IncludeOffset)
% Function for exponential fit.
% IncludeOffset: 0 neglect offset, 1 include offset

m = length(tData); % number of time-domain points
nExponents = length(x); % number of exponents

% Set up design matrix.
if IncludeOffset
  A = [exp(tData*x') ones(m,1)];
else
  A = exp(tData*x');
end

% Compute linear coefficients.
c = A\yData;

% Compute residuals.
Residuals = yData - A*c;

% Compute Jacobian
Jacobian = zeros(m,nExponents);
AA = A'*A;
for  j = 1:nExponents
  jj = -tData.*A(:,j);
  b = c(j)*(A'*jj);
  b(j) = b(j) - jj'*Residuals;
  Jacobian(:,j) = c(j)*jj - A*(AA\b);
end
   
