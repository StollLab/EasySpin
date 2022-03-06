% esfit_genetic   Genetic algorithm for least-squares fitting
%
%    x = esfit_genetic(fcn,nParams,Opt)
%
%    fcn      ... scalar function to minimize
%    nParams  ... number of parameters
%    Opt   ... options
%       .PopulationSize   number of individuals per generation
%       .EliteCount       number of elite individuals
%       .maxGenerations   maximum number of generations
%       .PrintLevel       1, if progress information should be printed
%       .TolFun           error threshold below which fitting stops

function bestx = esfit_genetic(fcn,lb,ub,Opt)

if nargin==0, help(mfilename); return; end

if nargin<3, Opt = struct; end

if ~isfield(Opt,'PopulationSize'), Opt.PopulationSize = 20; end
if ~isfield(Opt,'maxGenerations'), Opt.maxGenerations = 10000; end
if ~isfield(Opt,'EliteCount')
  Opt.EliteCount = max(2,ceil(0.1*Opt.PopulationSize));
end
if ~isfield(Opt,'PrintLevel'), Opt.PrintLevel = 0; end
if ~isfield(Opt,'Range'); Opt.Range = 1; end
if ~isfield(Opt,'TolFun'); Opt.TolFun = 1e-5; end
if ~isfield(Opt,'IterFcn'), Opt.IterFcn = []; end

stopCode = 0;

lb = lb(:).';
ub = ub(:).';
if numel(lb)~=numel(ub)
  error('Arrays for lower and upper bound must have the same number of elements.');
end
if any(lb>ub)
  error('Lower bounds must not be greater than upper bounds.');
end
nParams = numel(lb);

if Opt.PrintLevel
  fprintf('  %d parameters, range %g to %g\n',nParams,-Opt.Range,Opt.Range);
  fprintf('  population %d, elite %d\n',Opt.PopulationSize,Opt.EliteCount);
  fprintf('  %d generations\n',Opt.maxGenerations);
end

nEvals = 0;  % number of function evaluations

% Generate initial population
Population = lb+ (ub-lb).*rand(Opt.PopulationSize,nParams);

bestScore = inf;
bestx = zeros(size(Population(1,:)));

% Score initial population
if Opt.PrintLevel
  Opt.IterationPrintFunction('initial population');
end
Scores = ones(1,Opt.PopulationSize)*inf;
for k = 1:Opt.PopulationSize
  Scores(k) = fcn(Population(k,:));
  nEvals = nEvals+1;
  if Scores(k)<bestScore
    bestx = Population(k,:);
    bestScore = Scores(k);
  end
  info.bestx = bestx;
  info.minF = bestScore;
  info.nEvals = nEvals;
  info.iter = 0;
  if ~isempty(Opt.IterFcn)
    UserStop = Opt.IterFcn(info);
  else
    UserStop = false;
  end
end
[Scores,idx] = sort(Scores);
Population = Population(idx,:);
Fitness = Scores(end) - Scores;

gen = 1; % generation index

while true
  
  if stopCode, break; end

  if min(Scores)<bestScore, bestScore = min(Scores); end
  
  if Opt.PrintLevel
    str = sprintf('gen %5d:  min %g   mean %g',gen,min(Scores),mean(Scores));
    Opt.IterationPrintFunction(str);
  end
  
  if gen>=Opt.maxGenerations, stopCode = 1; break; end
  if bestScore<Opt.TolFun, stopCode = 2; break; end
  
  % (1) Selection
  %-----------------------------------------------
  RouletteWheel = [0 cumsum(Fitness)/sum(Fitness)];
  %Balls = rand(1,FitOpt.PopulationSize);
  Balls = (0.5+(0:Opt.PopulationSize-1))/Opt.PopulationSize;
  for k = 1:Opt.PopulationSize
    ParentsIdx(k) = sum(Balls(k)>RouletteWheel);
  end

  % (2) Recombination
  %-----------------------------------------------
  r = rand(2,Opt.PopulationSize);
  r = ParentsIdx(fix(Opt.PopulationSize*r+1));
  d = 0.1;
  a = rand(1,Opt.PopulationSize)*(1+2*d)-d;
  for k = 1:Opt.PopulationSize
    Offspring(k,:) = a(k)*Population(r(1,k),:) + (1-a(k))*Population(r(2,k));
  end
  
  % (3) Mutation
  %-----------------------------------------------
  InitialVariance = 0.3*(ub-lb);
  Variance = InitialVariance*(1-gen/Opt.maxGenerations);
  if Variance<0, Variance = 0; end
  Offspring = Offspring + randn(Opt.PopulationSize,nParams).*Variance;
  for p = 1:nParams
    Offspring(Offspring(:,p)<lb,p) = lb(p);
    Offspring(Offspring(:,p)>ub,p) = ub(p);
  end
  
  % (4) Reinsertion
  %-----------------------------------------------
  
  % Score offspring
  for k = 1:Opt.PopulationSize
    offScores(k) = fcn(Offspring(k,:));
    nEvals = nEvals+1;
    info.bestx = bestx;
    info.minF = bestScore;
    info.nEvals = nEvals;
    info.iter = gen;
    if ~isempty(Opt.IterFcn)
      UserStop = Opt.IterFcn(info);
    else
      UserStop = false;
    end
    if UserStop, stopCode = 3; break; end
  end
  if stopCode==3, break; end
  [offScores,idx] = sort(offScores);
  Offspring = Offspring(idx,:);

  idx = 1:(Opt.PopulationSize-Opt.EliteCount);
  Population(idx+Opt.EliteCount,:) = Offspring(idx,:);
  Scores(idx+Opt.EliteCount) = offScores(idx);

  [Scores,idx] = sort(Scores);
  Population = Population(idx,:);
  if Scores(1)<bestScore
    bestScore = Scores(1);
    bestx = Population(1,:);
  end
  Fitness = Scores(end) - Scores;

  info.bestx = bestx;
  info.minF = bestScore;
  info.nEvals = nEvals;
  info.iter = gen;
  if ~isempty(Opt.IterFcn)
    UserStop = Opt.IterFcn(info);
  else
    UserStop = false;
  end
  
  gen = gen + 1;
end

if Opt.PrintLevel>1
  switch stopCode
    case 1, msg = sprintf('Maximum number of generations (%d) reached.',Opt.maxGenerations);
    case 2, msg = sprintf('Error below threshold %g.',Opt.TolFun);
    case 3, msg = sprintf('Stopped by user.');
  end
  fprintf('Terminated: %s\n',msg);
end

return
