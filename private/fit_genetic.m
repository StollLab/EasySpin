% fit_genetic   Genetic algorithm for least-squares fitting
%
%    x = fit_genetic(funfcn,nParams,FitOpt,varargin)
%
%    funfcn ... scalar function to minimize
%    nParams ... number of parameters
%    FitOpt ... options
%       PopulationSize   number of individuals per generation
%       EliteCount       number of elite individuals
%       Range            parameter range (from -Range to +Range)
%       maxGenerations   maximum number of generations
%       PrintLevel          1, if progress information should be printed
%       TolFun           error threshold below which fitting stops

function bestx = fit_genetic(funfcn,nParams,FitOpt,varargin)

if (nargin==0), help(mfilename); return; end

global UserCommand;

if (nargin<3), FitOpt = struct('unused',NaN); end

if ~isfield(FitOpt,'PopulationSize'), FitOpt.PopulationSize = 20; end
if ~isfield(FitOpt,'maxGenerations'), FitOpt.maxGenerations = 40; end
if ~isfield(FitOpt,'EliteCount')
  FitOpt.EliteCount = max(2,ceil(0.1*FitOpt.PopulationSize));
end
if ~isfield(FitOpt,'PrintLevel'), FitOpt.PrintLevel = 0; end
if ~isfield(FitOpt,'Range'); FitOpt.Range = 1; end
if ~isfield(FitOpt,'TolFun'); FitOpt.TolFun = 0; end

if (FitOpt.PrintLevel)
  fprintf('  %d parameters, range %g to %g\n',nParams,-FitOpt.Range,FitOpt.Range);
  fprintf('  population %d, elite %d\n',FitOpt.PopulationSize,FitOpt.EliteCount);
  fprintf('  %d generations\n',FitOpt.maxGenerations);
end

stopCode = 0;

% Generate initial population
Population = FitOpt.Range*(2*rand(FitOpt.PopulationSize,nParams) - 1);

BestScore = inf;
bestx = Population(1,:)*0;


% Score initial population
Scores = ones(1,FitOpt.PopulationSize)*inf;
for k = 1:FitOpt.PopulationSize
  Scores(k) = feval(funfcn,Population(k,:),varargin{:});
  hLogLine = findobj('Tag','logLine');
  if Scores(k)<BestScore
    bestx = Population(k,:);
    BestScore = Scores(k);
  end
  if (UserCommand==1), stopCode = 3; break; end
  if (UserCommand==3)
    UserCommand = 0;
    set(FitOpt.hButtons(3),'Visible','off');
    set(FitOpt.hButtons(2),'Visible','on');
    if FitOpt.PrintLevel
      fprintf('       local:');
    end
    FitOpt2 = FitOpt;
    FitOpt2.PrintLevel = 0;
    [bestx,info] = FitOpt.fitlocal(funfcn,bestx,FitOpt2,varargin{1:end-1},FitOpt2);
    BestScore = info.F;
    if (UserCommand==4 || UserCommand==99), UserCommand = 0; end
    if FitOpt.PrintLevel
      fprintf('   error %g\n',info.F);
    end
    set(FitOpt.hButtons(2),'Visible','off');
    set(FitOpt.hButtons(3),'Visible','on');
  end
end
[Scores,idx] = sort(Scores);
Population = Population(idx,:);
Fitness = Scores(end) - Scores;

g = 1; % generation index

while 1
  
  if stopCode, break; end

  if (min(Scores)<BestScore), BestScore = min(Scores); end
  
  if FitOpt.PrintLevel
    str = sprintf('%5d:  min %g   mean %g',g,min(Scores),mean(Scores));
    if isempty(hLogLine)
      disp(str)
    else
      set(hLogLine,'String',str);
    end
  end
  
  if (g>=FitOpt.maxGenerations), stopCode = 1; break; end
  if (BestScore<FitOpt.TolFun), stopCode = 2; break; end
  if (UserCommand==1), stopCode = 3; break; end
  
  % (1) Selection
  %-----------------------------------------------
  RouletteWheel = [0 cumsum(Fitness)/sum(Fitness)];
  %Balls = rand(1,FitOpt.PopulationSize);
  Balls = (0.5+(0:FitOpt.PopulationSize-1))/FitOpt.PopulationSize;
  for k=1:FitOpt.PopulationSize
    ParentsIdx(k) = sum(Balls(k)>RouletteWheel);
  end

  % (2) Recombination
  %-----------------------------------------------
  r = rand(2,FitOpt.PopulationSize);
  r = ParentsIdx(fix(FitOpt.PopulationSize*r+1));
  d = 0.1;
  a = rand(1,FitOpt.PopulationSize)*(1+2*d)-d;
  for k=1:FitOpt.PopulationSize
    Offspring(k,:) = a(k)*Population(r(1,k),:) + (1-a(k))*Population(r(2,k));
  end
  
  % (3) Mutation
  %-----------------------------------------------
  InitialVariance = 0.3*2*FitOpt.Range;
  Variance = InitialVariance*(1-g/FitOpt.maxGenerations);
  if (Variance<0), Variance = 0; end
  for k=1:FitOpt.PopulationSize
    Offspring(k,:) = Offspring(k,:) + randn(1,nParams)*Variance;
  end
  Offspring(Offspring<-FitOpt.Range) = -FitOpt.Range;
  Offspring(Offspring>+FitOpt.Range) = +FitOpt.Range;
  
  % (4) Reinsertion
  %-----------------------------------------------

  % Score offspring
  for k = 1:FitOpt.PopulationSize
    offScores(k) = feval(funfcn,Offspring(k,:),varargin{:});
    if (UserCommand==1), stopCode = 3; break; end
    if (UserCommand==3)
      UserCommand = 0;
      set(FitOpt.hButtons(3),'Visible','off');
      set(FitOpt.hButtons(2),'Visible','on');
      if FitOpt.PrintLevel
        fprintf('       local:');
      end
      FitOpt2 = FitOpt;
      FitOpt2.PrintLevel = 0;
      [bestx,info] = fit_simplex(funfcn,bestx,FitOpt2,varargin{1:end-1},FitOpt2);
      BestScore = info.F;
      if (UserCommand==4), UserCommand = 0; end
      if FitOpt.PrintLevel
        fprintf('   error %g\n',info.F);
      end
      set(FitOpt.hButtons(2),'Visible','off');
      set(FitOpt.hButtons(3),'Visible','on');
    end
  end
  if (stopCode==3), break; end
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

if FitOpt.PrintLevel
  switch (stopCode)
    case 1, msg = sprintf('Maximum number of generations (%d) reached.',FitOpt.maxGenerations);
    case 2, msg = sprintf('Error below threshold %g.',FitOpt.TolFun);
    case 3, msg = sprintf('Stopped by user.');
  end
  fprintf('Terminated: %s\n',msg);
end

return
