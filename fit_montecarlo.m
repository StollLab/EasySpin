function [bestx,info] = fit_montecarlo(funfcn,nParameters,FitOpt,varargin)

if ~isfield(FitOpt,'N'); FitOpt.N = 10000; end
if ~isfield(FitOpt,'PrintLevel'); FitOpt.PrintLevel = 1; end
if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-5; end

nMonteCarlo = FitOpt.N;

if FitOpt.PrintLevel
  fprintf('Doing %d random trials...\n',nMonteCarlo);
end

global UserCommand

minerror = inf;
bestx = zeros(nParameters,1);
startTime = cputime;
  
stopCode = 0;
for k = 1:nMonteCarlo
  X = 2*rand(nParameters,1) - 1;
  thiserror = feval(funfcn,X,varargin{:});
  hLogLine = findobj('Tag','logLine');
  if (thiserror<minerror)
    minerror = thiserror;
    bestx = X;
    if FitOpt.PrintLevel
      str=sprintf('  Trial %4d:   error %0.5e  best so far',k,thiserror);
      if isempty(hLogLine)
        disp(str);
      else
        set(hLogLine,'String',str);
      end
    end
  end
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
    if (UserCommand==4) || (UserCommand==99), UserCommand = 0; end
    if FitOpt.PrintLevel
      fprintf('   error %0.5e\n',info.F);
    end
    set(FitOpt.hButtons(2),'Visible','off');
    set(FitOpt.hButtons(3),'Visible','on');
  end
  elapsedTime = (cputime-startTime)/60;
  if (elapsedTime>FitOpt.maxTime), stopCode = 1; break; end
  if (UserCommand==1), stopCode = 2; break; end
  if (thiserror<FitOpt.TolFun), stopCode = 3; break; end
end

if FitOpt.PrintLevel
  switch (stopCode)
    case 0, msg = 'Terminated: maximum number of trials reached.';
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = 'Terminated: Stopped by user.';
    case 3, msg = sprintf('Terminated: Found a parameter set with error less than %g.',FitOpt.TolFun);
  end
  disp(msg);
end
