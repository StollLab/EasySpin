function [gX,info] = esfit_particleswarm(funfcn,nParameters,FitOpt,varargin)

if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'PrintLevel'); FitOpt.PrintLevel = 1; end
if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-5; end
if ~isfield(FitOpt,'IterationPrintFunction') || ...
  isempty(FitOpt.IterationPrintFunction), FitOpt.IterationPrintFunction = @(str)0; end
if ~isfield(FitOpt,'nParticles'); FitOpt.nParticles = 30; end
if ~isfield(FitOpt,'SwarmParams'), FitOpt.SwarmParams = [0.2 0.5 2 1]; end

global UserCommand
if isempty(UserCommand), UserCommand = NaN; end

nParticles = FitOpt.nParticles;

if FitOpt.PrintLevel
  fprintf('Particle swarm size %d.\n',nParticles);
end

k = FitOpt.SwarmParams(1); % velocity clamping
w = FitOpt.SwarmParams(2); % inertial coefficient
c1 = FitOpt.SwarmParams(3); % cognitive coefficient
c2 = FitOpt.SwarmParams(4); % social coefficient

if FitOpt.PrintLevel
  fprintf('Particle swarm parameters: k = %g, w = %g, c1 = %g, c2 = %g.\n',k,w,c1,c2);
end

X = 2*rand(nParameters,nParticles) - 1;
bestX = X;
gX = bestX(:,1);
v = (2*rand(nParameters,nParticles)-1)*k;
besterror = ones(1,nParticles)*inf;
globalbesterror = inf;
minerror = inf;
startTime = cputime;
  
FitOpt.IterationPrintFunction('initial iteration');

iIteration = 1;
stopCode = 0;
while 1
  iIteration = iIteration + 1;
  
  for p = 1:nParticles
    if (UserCommand==1), stopCode = 2; break; end
    thiserror(p) = funfcn(X(:,p),varargin{:});
    if thiserror(p)<besterror(p)
      besterror(p) = thiserror(p);
      bestX(:,p) = X(:,p);
    end
  end
  
  [val,idx] = min(thiserror);
  if thiserror(idx)<globalbesterror
    globalbesterror = thiserror(idx);
    gX = bestX(:,idx);
  end
  
  for p = 1:nParticles
    v(:,p) = w*v(:,p) + c1*rand*(bestX(:,p)-X(:,p)) + c2*rand*(gX-X(:,p));
    X(:,p) = X(:,p) + v(:,p);
    X(:,p) = min(max(X(:,p),-1),+1);
  end
  
  if (globalbesterror<minerror)
    minerror = thiserror;
  end
  
  if FitOpt.PrintLevel
    str=sprintf('  Iteration %4d:   error %0.5e  best so far',iIteration,globalbesterror);
    FitOpt.IterationPrintFunction(str);
  end
  
  elapsedTime = (cputime-startTime)/60;
  if (elapsedTime>FitOpt.maxTime), stopCode = 1; break; end
  if (UserCommand==1), stopCode = 2; break; end
  if (globalbesterror<FitOpt.TolFun), stopCode = 3; break; end
end

if FitOpt.PrintLevel>1
  switch (stopCode)
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = 'Terminated: Stopped by user.';
    case 3, msg = sprintf('Terminated: Found a parameter set with error less than %g.',FitOpt.TolFun);
  end
  disp(msg);
end
