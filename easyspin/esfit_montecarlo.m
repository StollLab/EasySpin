function [bestx,info] = esfit_montecarlo(funfcn,nParameters,FitOpt,varargin)

if ~isfield(FitOpt,'nTrials'); FitOpt.nTrials = 20000; end
if ~isfield(FitOpt,'PrintLevel'); FitOpt.PrintLevel = 1; end
if ~isfield(FitOpt,'maxTime'), FitOpt.maxTime = inf; end
if ~isfield(FitOpt,'TolFun'), FitOpt.TolFun = 1e-5; end
if ~isfield(FitOpt,'IterationPrintFunction') || ...
  isempty(FitOpt.IterationPrintFunction)
  FitOpt.IterationPrintFunction = @iterationprint;
end

global UserCommand
if isempty(UserCommand), UserCommand = NaN; end

minerror = inf;
bestx = zeros(nParameters,1);
startTime = cputime;
  
stopCode = 0;
for k = 1:FitOpt.nTrials
  
  X = 2*rand(nParameters,1) - 1;
  thiserror = feval(funfcn,X,varargin{:});
  
  if (thiserror<minerror)
    minerror = thiserror;
    bestx = X;
    if FitOpt.PrintLevel
      str = sprintf('%d:   error %0.5e  best so far',k,thiserror);
      FitOpt.IterationPrintFunction(str);
    end
  end
    
  elapsedTime = (cputime-startTime)/60;
  if (elapsedTime>FitOpt.maxTime), stopCode = 1; break; end
  if (UserCommand==1), stopCode = 2; break; end
  if (thiserror<FitOpt.TolFun), stopCode = 3; break; end

end

if FitOpt.PrintLevel>1
  switch (stopCode)
    case 0, msg = 'Terminated: maximum number of trials reached.';
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = 'Terminated: Stopped by user.';
    case 3, msg = sprintf('Terminated: Found a parameter set with error less than %g.',FitOpt.TolFun);
  end
  disp(msg);
end

end

function iterationprint(str)
hLogLine = findobj('Tag','logLine');
if isempty(hLogLine)
  disp(str);
else
  set(hLogLine,'String',str);
end
end