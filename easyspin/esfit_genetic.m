% esfit_genetic   Genetic algorithm for least-squares fitting
%
%    x = esfit_genetic(fcn,nParams,FitOpt)
%
%    fcn      ... scalar function to minimize
%    nParams  ... number of parameters
%    FitOpt   ... options
%       .PopulationSize   number of individuals per generation
%       .EliteCount       number of elite individuals
%       .maxGenerations   maximum number of generations
%       .Verbosity       1, if progress information should be printed
%       .TolFun           error threshold below which fitting stops

function [bestx,info] = esfit_genetic(fcn,lb,ub,FitOpt)

if nargin==0, help(mfilename); return; end
if nargin<3
  error('At least 3 inputs expected (function, lb, ub).');
end
if nargin==3, FitOpt = struct; end

DefOpt = esfit_algdefaults('genetic algorithm');
FitOpt = adddefaults(FitOpt,DefOpt);

if ~isfield(FitOpt,'EliteCount') || isempty(FitOpt.EliteCount)
  FitOpt.EliteCount = max(2,ceil(0.1*FitOpt.PopulationSize));
end

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

if FitOpt.Verbosity
  msg{1} = sprintf('  %d parameters',nParams);
  msg{2} = sprintf('  population %d, elite %d',FitOpt.PopulationSize,FitOpt.EliteCount);
  msg{3} = sprintf('  %d generations',FitOpt.maxGenerations);
  FitOpt.InfoPrintFunction(msg);
end

nEvals = 0;  % number of function evaluations
startTime = cputime;

% Generate initial population
Population = lb+ (ub-lb).*rand(FitOpt.PopulationSize,nParams);

bestScore = inf;
bestx = zeros(size(Population(1,:)));

% Score initial population
if FitOpt.Verbosity
  FitOpt.IterationPrintFunction('initial population');
end
Scores = ones(1,FitOpt.PopulationSize)*inf;
for k = 1:FitOpt.PopulationSize
  Scores(k) = fcn(Population(k,:));
  nEvals = nEvals+1;
  newbest = Scores(k)<bestScore;
  if newbest
    bestx = Population(k,:);
    bestScore = Scores(k);
  end
  info.bestx = bestx;
  info.minF = bestScore;
  info.nEvals = nEvals;
  info.iter = 0;
  info.newbest = newbest;
  if ~isempty(FitOpt.IterFcn)
    stopCode = FitOpt.IterFcn(info);
    if stopCode, break; end
  else
    stopCode = false;
  end
end
[Scores,idx] = sort(Scores);
Population = Population(idx,:);
Fitness = Scores(end) - Scores;

gen = 1; % generation index

while true
  
  if stopCode, break; end

  if min(Scores)<bestScore, bestScore = min(Scores); end
  
  if FitOpt.Verbosity
    str = sprintf('gen %5d:  min %g   mean %g',gen,min(Scores),mean(Scores));
    FitOpt.IterationPrintFunction(str);
  end
  
  if gen>=FitOpt.maxGenerations, stopCode = 0; break; end
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>FitOpt.maxTime, stopCode = 1; break; end
  if bestScore<FitOpt.TolFun, stopCode = 2; break; end
  
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
    Offspring(k,:) = a(k)*Population(r(1,k),:) + (1-a(k))*Population(r(2,k),:);
  end
  
  % (3) Mutation
  %-----------------------------------------------
  % Mutation with consideration of limits
  % (to avoid excessive weight of lb and ub values)
  % based on: Blume C., Wilfried J.
  % "GLEAM - An Evolutionary Algorithm for Planning 
  % and Control Based on Evolution Strategy." 
  % Conf. Proc. of Genetic and Evolutionary Computation Conference
  % Late Breaking Papers (2002)
  nk = 10; % number of sections
  r = 2*randi([0 1],1,FitOpt.PopulationSize)-1;
  j = randi(nk-1,1,FitOpt.PopulationSize); % select random section
  for k = 1:FitOpt.PopulationSize
    if r(k)==-1 % decrease
      rangelim = [lb; Offspring(k,:)];
    elseif r(k)==1 % increase
      rangelim = [Offspring(k,:); ub];
    end
    delta = (rangelim(2,:)-rangelim(1,:))/nk;
    Offspring(k,:) = Offspring(k,:) + r(k)*rand(1,nParams).*delta*j(k);
  end
  for p = 1:nParams
    Offspring(Offspring(:,p)<lb(p),p) = lb(p);
    Offspring(Offspring(:,p)>ub(p),p) = ub(p);
  end
  
  % (4) Reinsertion
  %-----------------------------------------------
  
  % Score offspring
  for k = 1:FitOpt.PopulationSize
    offScores(k) = fcn(Offspring(k,:));
    nEvals = nEvals+1;
    info.nEvals = nEvals;
    info.iter = gen;
    if ~isempty(FitOpt.IterFcn)
      UserStop = FitOpt.IterFcn(info);
    else
      UserStop = false;
    end
    if UserStop, stopCode = 3; break; end
  end
  if stopCode==3, break; end
  [offScores,idx] = sort(offScores);
  Offspring = Offspring(idx,:);

  idx = 1:(FitOpt.PopulationSize-FitOpt.EliteCount);
  Population(idx+FitOpt.EliteCount,:) = Offspring(idx,:);
  Scores(idx+FitOpt.EliteCount) = offScores(idx);

  [Scores,idx] = sort(Scores);
  Population = Population(idx,:);
  newbest = Scores(1)<bestScore;
  if newbest
    bestScore = Scores(1);
    bestx = Population(1,:);
  end
  Fitness = Scores(end) - Scores;

  info.bestx = bestx;
  info.minF = bestScore;
  info.nEvals = nEvals;
  info.iter = gen;
  info.newbest = newbest;
  if ~isempty(FitOpt.IterFcn)
    UserStop = FitOpt.IterFcn(info);
  else
    UserStop = false;
  end
  if UserStop, stopCode = 3; break; end
  
  gen = gen + 1;
end

if FitOpt.Verbosity>0
  switch stopCode
    case 0, msg = sprintf('Maximum number of generations (%d) reached.',FitOpt.maxGenerations);
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = sprintf('Error below threshold %g.',FitOpt.TolFun);
    case 3, msg = sprintf('Stopped by user.');
  end
  FitOpt.InfoPrintFunction(sprintf('Terminated: %s\n',msg));
end

return
