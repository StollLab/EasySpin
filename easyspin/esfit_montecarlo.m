function [bestx,info] = esfit_montecarlo(funfcn,nParams,FitOpt,varargin)

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

lb = -ones(nParams,1);
ub = +ones(nParams,1);
if any(lb>ub)
  error('Lower bounds must not be greater than upper bounds.');
end

Fmin = inf;
bestx = zeros(nParams,1);
startTime = cputime;
  
stopCode = 0;
for iTrial = 1:FitOpt.nTrials
  
  X = lb + (ub-lb).*rand(nParams,1);
  F = feval(funfcn,X,varargin{:});
  
  if F<Fmin
    Fmin = F;
    bestx = X;
    if FitOpt.PrintLevel
      str = sprintf('%d:   error %0.5e  best so far',iTrial,F);
      FitOpt.IterationPrintFunction(str);
    end
  end
    
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>FitOpt.maxTime, stopCode = 1; end
  if UserCommand==1, stopCode = 2; end
  if F<FitOpt.TolFun, stopCode = 3; end
  
  if stopCode, break; end
  
end

if FitOpt.PrintLevel>1
  switch stopCode
    case 0, msg = 'Terminated: maximum number of trials reached.';
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = 'Terminated: Stopped by user.';
    case 3, msg = sprintf('Terminated: Found a parameter set with error less than %g.',FitOpt.TolFun);
  end
  disp(msg);
end

info.nTrials = iTrial;
info.elapsedTime = elapsedTime;

end

function iterationprint(str)
hLogLine = findobj('Tag','logLine');
if isempty(hLogLine)
  disp(str);
else
  set(hLogLine,'String',str);
end
end
