% esfit_genetic   Genetic algorithm for least-squares fitting
%
%    x = esfit_genetic(fcn,nParams,FitOpt,varargin)
%
%    fcn      ... scalar function to minimize
%    nParams  ... number of parameters
%    FitOpt   ... options
%       .PopulationSize   number of individuals per generation
%       .EliteCount       number of elite individuals
%       .maxGenerations   maximum number of generations
%       .PrintLevel       1, if progress information should be printed
%       .TolFun           error threshold below which fitting stops

function bestx = esfit_genetic(fcn,lb,ub,FitOpt)

if nargin==0, help(mfilename); return; end

global UserCommand
if isempty(UserCommand), UserCommand = NaN; end

if nargin<3, FitOpt = struct; end

if ~isfield(FitOpt,'PopulationSize'), FitOpt.PopulationSize = 20; end
if ~isfield(FitOpt,'maxGenerations'), FitOpt.maxGenerations = 10000; end
if ~isfield(FitOpt,'EliteCount')
  FitOpt.EliteCount = max(2,ceil(0.1*FitOpt.PopulationSize));
end
if ~isfield(FitOpt,'PrintLevel'), FitOpt.PrintLevel = 0; end
if ~isfield(FitOpt,'Range'); FitOpt.Range = 1; end
if ~isfield(FitOpt,'TolFun'); FitOpt.TolFun = 1e-5; end

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

if FitOpt.PrintLevel
  fprintf('  %d parameters, range %g to %g\n',nParams,-FitOpt.Range,FitOpt.Range);
  fprintf('  population %d, elite %d\n',FitOpt.PopulationSize,FitOpt.EliteCount);
  fprintf('  %d generations\n',FitOpt.maxGenerations);
end

% Generate initial population
Population = lb+ (ub-lb).*rand(FitOpt.PopulationSize,nParams);

BestScore = inf;
bestx = zeros(size(Population(1,:)));

% Score initial population
if FitOpt.PrintLevel
  FitOpt.IterationPrintFunction('initial population');
end
Scores = ones(1,FitOpt.PopulationSize)*inf;
for k = 1:FitOpt.PopulationSize
  Scores(k) = feval(fcn,Population(k,:));
  if Scores(k)<BestScore
    bestx = Population(k,:);
    BestScore = Scores(k);
  end
  if UserCommand==1, stopCode = 3; break; end
end
[Scores,idx] = sort(Scores);
Population = Population(idx,:);
Fitness = Scores(end) - Scores;

g = 1; % generation index

while true
  
  if stopCode, break; end

  if min(Scores)<BestScore, BestScore = min(Scores); end
  
  if FitOpt.PrintLevel
    str = sprintf('gen %5d:  min %g   mean %g',g,min(Scores),mean(Scores));
    FitOpt.IterationPrintFunction(str);
  end
  
  if g>=FitOpt.maxGenerations, stopCode = 1; break; end
  if BestScore<FitOpt.TolFun, stopCode = 2; break; end
  if UserCommand==1, stopCode = 3; break; end
  
  % (1) Selection
  %-----------------------------------------------
  RouletteWheel = [0 cumsum(Fitness)/sum(Fitness)];
  %Balls = rand(1,FitOpt.PopulationSize);
  Balls = (0.5+(0:FitOpt.PopulationSize-1))/FitOpt.PopulationSize;
  for k = 1:FitOpt.PopulationSize
    ParentsIdx(k) = sum(Balls(k)>RouletteWheel);
  end

  % (2) Recombination
  %-----------------------------------------------
  r = rand(2,FitOpt.PopulationSize);
  r = ParentsIdx(fix(FitOpt.PopulationSize*r+1));
  d = 0.1;
  a = rand(1,FitOpt.PopulationSize)*(1+2*d)-d;
  for k = 1:FitOpt.PopulationSize
    Offspring(k,:) = a(k)*Population(r(1,k),:) + (1-a(k))*Population(r(2,k));
  end
  
  % (3) Mutation
  %-----------------------------------------------
  InitialVariance = 0.3*(ub-lb);
  Variance = InitialVariance*(1-g/FitOpt.maxGenerations);
  if Variance<0, Variance = 0; end
  for k=1:FitOpt.PopulationSize
    Offspring(k,:) = Offspring(k,:) + randn(1,nParams)*Variance;
  end
  Offspring(Offspring<lb) = lb;
  Offspring(Offspring>ub) = ub;
  
  % (4) Reinsertion
  %-----------------------------------------------
  
  % Score offspring
  for k = 1:FitOpt.PopulationSize
    offScores(k) = feval(fcn,Offspring(k,:));
    if UserCommand==1, stopCode = 3; break; end
  end
  if stopCode==3, break; end
  [offScores,idx] = sort(offScores);
  Offspring = Offspring(idx,:);

  idx = 1:(FitOpt.PopulationSize-FitOpt.EliteCount);
  Population(idx+FitOpt.EliteCount,:) = Offspring(idx,:);
  Scores(idx+FitOpt.EliteCount) = offScores(idx);

  [Scores,idx] = sort(Scores);
  Population = Population(idx,:);
  if Scores(1)<BestScore
    BestScore = Scores(1);
    bestx = Population(1,:);
  end
  Fitness = Scores(end) - Scores;
  
  g = g + 1;
end

if FitOpt.PrintLevel>1
  switch stopCode
    case 1, msg = sprintf('Maximum number of generations (%d) reached.',FitOpt.maxGenerations);
    case 2, msg = sprintf('Error below threshold %g.',FitOpt.TolFun);
    case 3, msg = sprintf('Stopped by user.');
  end
  fprintf('Terminated: %s\n',msg);
end

return
