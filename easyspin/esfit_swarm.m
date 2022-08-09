function [gX,info] = esfit_swarm(fcn,lb,ub,FitOpt)

if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'Verbosity'); FitOpt.Verbosity = 1; end
if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-5; end
if ~isfield(FitOpt,'IterationPrintFunction'), FitOpt.IterationPrintFunction = @(str)str; end
if ~isfield(FitOpt,'InfoPrintFunction'), FitOpt.InfoPrintFunction = @(str)str; end
if ~isfield(FitOpt,'IterFcn'), FitOpt.IterFcn = []; end
if ~isfield(FitOpt,'SwarmParams'), FitOpt.SwarmParams = [0.2 0.5 2 1]; end
if ~isfield(FitOpt,'TolStallIter'), FitOpt.TolStallIter = 6; end

global UserCommand
if isempty(UserCommand), UserCommand = NaN; end

lb = lb(:);
ub = ub(:);
if numel(lb)~=numel(ub)
  error('Arrays for lower and upper bound must have the same number of elements.');
end
if any(lb>ub)
  error('Lower bounds must not be greater than upper bounds.');
end
nParams = numel(lb);

if ~isfield(FitOpt,'nParticles'); FitOpt.nParticles = 20 + nParams*10; end

nParticles = FitOpt.nParticles;

k = FitOpt.SwarmParams(1); % velocity clamping
w = FitOpt.SwarmParams(2); % inertial coefficient
c1 = FitOpt.SwarmParams(3); % cognitive coefficient
c2 = FitOpt.SwarmParams(4); % social coefficient

if FitOpt.Verbosity
  msg{1} = sprintf('Particle swarm optimization parameters:');
  msg{2} = sprintf('   n = %d (number of particles)',nParticles);
  msg{3} = sprintf('   k = %g (velocity clampling)',k);
  msg{4} = sprintf('   k = %g (inertial coefficient)',w);
  msg{5} = sprintf('   c1 = %g (cognitive coefficient)',c1);
  msg{6} = sprintf('   c2 = %g (social coefficient)',c2);
  FitOpt.InfoPrintFunction(msg);
end

if FitOpt.Verbosity
  FitOpt.InfoPrintFunction(sprintf('Initializing swarm...'));
end

% Initial particle positions and velocities
X = lb + (ub-lb).*rand(nParams,nParticles);
v = (ub-lb).*rand(nParams,nParticles)*k;

bestX = X;
gX = bestX(:,1);
bestF = inf(1,nParticles);
F = inf(1,nParticles);
globalbestF = inf;
startTime = cputime;
nStalledIterations = 0; % counts the number of iterations globalbestF hasn't changed

if FitOpt.Verbosity
  FitOpt.InfoPrintFunction(sprintf('Starting iterations...'));
end
iIteration = 0;
stopCode = 0;
while stopCode==0
  
  iIteration = iIteration + 1;
  
  % Evaluate functions for all particles
  for p = 1:nParticles
    if UserCommand==1, stopCode = 2; break; end
    F(p) = fcn(X(:,p));
  end
  
  % Update best fit so far for each particle
  idx = F<bestF;
  bestF(idx) = F(idx);
  bestX(:,idx) = X(:,idx);
  
  % Find overall best fit so far
  [~,idx] = min(F);
  newbest = F(idx)<globalbestF;
  if newbest
    nStalledIterations = 0;
    globalbestF = F(idx);
    gX = X(:,idx);
  else
    nStalledIterations = nStalledIterations+1;
  end
  
  % Move all particles
  for p = 1:nParticles
    v(:,p) = w*v(:,p) + c1*rand*(bestX(:,p)-X(:,p)) + c2*rand*(gX-X(:,p));
    X(:,p) = X(:,p) + v(:,p);
    X(:,p) = min(max(X(:,p),lb),ub); % constrain to [lb ub]
  end
  
  if FitOpt.Verbosity
    str = sprintf('  Iteration %4d:   value %0.5e  best so far (%d)',iIteration,globalbestF,nStalledIterations);
    FitOpt.IterationPrintFunction(str);
  end
  
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>FitOpt.maxTime, stopCode = 1; end
  if UserCommand==1, stopCode = 2; end
  if globalbestF<FitOpt.TolFun, stopCode = 3; end
  if nStalledIterations==FitOpt.TolStallIter, stopCode = 4; end
  if ~isempty(FitOpt.IterFcn)
    info.newbest = newbest;
    UserStop = FitOpt.IterFcn(info);
    if UserStop
      stopCode = 2;
    end
  end
  
end

if FitOpt.Verbosity>1
  switch stopCode
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = 'Terminated: Stopped by user.';
    case 3, msg = sprintf('Terminated: Found a parameter set with function value less than %g.',FitOpt.TolFun);
    case 4, msg = sprintf('Terminated: Function value didn''t change for %d iterations.',FitOpt.TolStallIter);
  end
  FitOpt.InfoPrintFunction(msg);
end

info.nIterations = iIteration;
info.elapsedTime = elapsedTime;
info.newbest = newbest;
end
