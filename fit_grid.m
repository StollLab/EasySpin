function [bestx,info] = fit_grid(funfcn,nParameters,FitOpt,varargin)

global UserCommand

GridSize = FitOpt.GridSize;
if numel(GridSize)==1
  GridSize = GridSize*ones(1,nParameters);
end
if numel(GridSize)~=nParameters
  error('FitOpt.GridSize must have as many elements as there are fitting parameters.');
end
if any(GridSize<1)
  error('At least one grid point per parameter is needed.');
end

nGridPoints = prod(GridSize);
if nGridPoints>FitOpt.maxGridPoints
  error('Cannot do grid search with more than %d points. Reduce number of parameters.',FitOpt.maxGridPoints);
end

for p = 1:nParameters
  if GridSize(p)==1
    grid{p} = 0;
  else
    grid{p} = linspace(-1,1,GridSize(p));
  end
end

X = cell(1,nParameters);
[X{:}] = ndgrid(grid{:});
for k=1:nParameters, X{k} = X{k}(:); end
X = [X{end:-1:1}];

FitOpt.RandomizeGrid = 1;
if FitOpt.RandomizeGrid, X = X(randperm(nGridPoints),:); end

minerror = inf;
bestx = zeros(nParameters,1);
startTime = cputime;

if FitOpt.PrintLevel
  fprintf('%d parameters, %d grid points total\n',...
    nParameters,nGridPoints);
end

stopCode = 0;
for k = 1:nGridPoints
  thiserror = feval(funfcn,X(k,:),varargin{:});
  hLogLine = findobj('Tag','logLine');
  if (thiserror<minerror)
    minerror = thiserror;
    bestx = X(k,:);
    if FitOpt.PrintLevel
      str = sprintf('  Point %4d:   error %0.5e  best so far',k,thiserror);
      if isempty(hLogLine)
        disp(str)
      else
        set(hLogLine,'String',str);
      end
    end
  end
  if FitOpt.Plot
    if (UserCommand==3)
      set(FitOpt.hButtons(3),'Visible','off');
      set(FitOpt.hButtons(2),'Visible','on');
      UserCommand = 0;
      if FitOpt.PrintLevel
        fprintf('       local:');
      end
      FitOpt2 = FitOpt;
      FitOpt2.PrintLevel = 0;
      FitOpt2.TolX = FitOpt2.TolStep;
      [bestx,info] = FitOpt.fitlocal(funfcn,bestx,FitOpt2,varargin{1:end-1},FitOpt2);
      if (UserCommand==4 || UserCommand == 99), UserCommand = 0; end
      if FitOpt.PrintLevel
        fprintf('   error %0.5e\n',info.F);
      end
      set(FitOpt.hButtons(2),'Visible','off');
      set(FitOpt.hButtons(3),'Visible','on');
    end
  end
  elapsedTime = (cputime-startTime)/60;
  if (elapsedTime>FitOpt.maxTime), stopCode = 1; break; end
  if (UserCommand==1), stopCode = 2; break; end
  if (thiserror<FitOpt.TolFun), stopCode = 3; break; end

end

switch (stopCode)
  case 0, msg = 'Terminated: all grid points searched.';
  case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
  case 2, msg = 'Terminated: Stopped by user.';
  case 3, msg = sprintf('Terminated: Found a parameter set with error less than %g.',FitOpt.TolFun);
end

if FitOpt.PrintLevel, disp(msg); end
