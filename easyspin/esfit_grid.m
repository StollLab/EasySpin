function [bestx,info] = esfit_grid(fcn,lb,ub,FitOpt)

global UserCommand
if isempty(UserCommand), UserCommand = NaN; end

if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-5; end
if ~isfield(FitOpt,'GridSize'), FitOpt.GridSize = 7; end
if ~isfield(FitOpt,'RandomizeGrid'), FitOpt.RandomizeGrid = true; end

lb = lb(:);
ub = ub(:);
if numel(lb)~=numel(ub)
  error('Arrays for lower and upper bound must have the same number of elements.');
end
if any(lb>ub)
  error('Lower bounds must not be greater than upper bounds.');
end
nParams = numel(lb);

GridSize = FitOpt.GridSize;
if numel(GridSize)==1
  GridSize = GridSize*ones(1,nParams);
end
if numel(GridSize)~=nParams
  error('FitOpt.GridSize must have as many elements as there are fitting parameters.');
end
if any(GridSize<2)
  error('At least one grid point per parameter is needed.');
end

nGridPoints = prod(GridSize);
if nGridPoints>FitOpt.maxGridPoints
  error('Cannot do grid search with more than %d points. Reduce number of parameters.',FitOpt.maxGridPoints);
end

for p = 1:nParams
  grid{p} = linspace(lb(p),ub(p),GridSize(p));
end

X = cell(1,nParams);
[X{:}] = ndgrid(grid{:});
for k=1:nParams, X{k} = X{k}(:); end
X = [X{:}];

if FitOpt.RandomizeGrid
  X = X(randperm(nGridPoints),:);
end

minF = inf;
bestx = zeros(nParams,1);
startTime = cputime;

if FitOpt.PrintLevel
  fprintf('%d parameters, %d grid points total\n',...
    nParams,nGridPoints);
end

stopCode = 0;
for k = 1:nGridPoints
  
  F = fcn(X(k,:));
  
  if F<minF
    minF = F;
    bestx = X(k,:);
    if FitOpt.PrintLevel
      str = sprintf('  Point %4d/%d:   error %0.5e  best so far',k,nGridPoints,F);
      FitOpt.IterationPrintFunction(str);
    end
  end
  
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>FitOpt.maxTime, stopCode = 1; end
  if UserCommand==1, stopCode = 2; end
  if F<FitOpt.TolFun, stopCode = 3; end
  
  if stopCode, break; end
  
end

switch stopCode
  case 0, msg = 'Terminated: all grid points searched.';
  case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
  case 2, msg = 'Terminated: Stopped by user.';
  case 3, msg = sprintf('Terminated: Found a parameter set with error less than %g.',FitOpt.TolFun);
end

if FitOpt.PrintLevel>1, disp(msg); end

return
