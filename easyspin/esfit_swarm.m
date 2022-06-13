function [gX,info] = esfit_swarm(fcn,lb,ub,FitOpt)

if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'Verbosity'); FitOpt.Verbosity = 1; end
if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-5; end
if ~isfield(FitOpt,'IterationPrintFunction'), FitOpt.IterationPrintFunction = []; end
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
  fprintf('Particle swarm optimization parameters:\n');
  fprintf('   n = %d (number of particles)\n',nParticles);
  fprintf('   k = %g (velocity clampling)\n',k);
  fprintf('   k = %g (inertial coefficient)\n',w);
  fprintf('   c1 = %g (cognitive coefficient)\n',c1);
  fprintf('   c2 = %g (social coefficient)\n',c2);
end

if FitOpt.Verbosity
  fprintf('Initializing swarm...\n');
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
  fprintf('Starting iterations...\n');
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
  if F(idx)<globalbestF
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
  
  if FitOpt.Verbosity && ~isempty(FitOpt.IterationPrintFunction)
    str = sprintf('  Iteration %4d:   value %0.5e  best so far (%d)',iIteration,globalbestF,nStalledIterations);
    FitOpt.IterationPrintFunction(str);
  end
  
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>FitOpt.maxTime, stopCode = 1; end
  if UserCommand==1, stopCode = 2; end
  if globalbestF<FitOpt.TolFun, stopCode = 3; end
  if nStalledIterations==FitOpt.TolStallIter, stopCode = 4; end
  
end

if FitOpt.Verbosity>1
  switch stopCode
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = 'Terminated: Stopped by user.';
    case 3, msg = sprintf('Terminated: Found a parameter set with function value less than %g.',FitOpt.TolFun);
    case 4, msg = sprintf('Terminated: Function value didn''t change for %d iterations.',FitOpt.TolStallIter);
  end
  disp(msg);
end

info.nIterations = iIteration;
info.elapsedTime = elapsedTime;

end
