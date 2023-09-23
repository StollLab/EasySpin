% esfit_montecarlo   Monte Carlo algorithm for least-squares fitting
%
%    x = esfit_montecarlo(fcn,lb,ub,FitOpt)
%
%    fcn ...  scalar function to minimize, f(x), where x is an array
%    lb   ... lower bounds for x
%    ub   ... upper bounds for x
%    FitOpt ... options
%       nTrials       number of trials
%       maxTime       maximum time to run, in minutes
%       Verbosity    1, if progress information should be printed
%       TolFun        error threshold below which fitting stops

function [bestx,info] = esfit_montecarlo(fcn,lb,ub,FitOpt)

if nargin==0, help(mfilename); return; end
if nargin<3
  error('At least 3 inputs expected (function, lb, ub).');
end
if nargin==3, FitOpt = struct; end

% Set options
DefOpt = esfit_algdefaults('Monte Carlo');
FitOpt = adddefaults(FitOpt,DefOpt);

lb = lb(:);
ub = ub(:);
if numel(lb)~=numel(ub)
  error('Arrays for lower and upper bound must have the same number of elements.');
end
if any(lb>ub)
  error('Lower bounds must not be greater than upper bounds.');
end
nParams = numel(lb);

Fmin = inf;
bestx = zeros(nParams,1);
startTime = cputime;
  
stopCode = 0;
for iTrial = 1:FitOpt.nTrials
  
  X = lb + (ub-lb).*rand(nParams,1);
  F = fcn(X);
  
  newbest = F<Fmin;
  if newbest
    Fmin = F;
    bestx = X;
    if FitOpt.Verbosity
      str = sprintf('%d:   error %0.5e  best so far',iTrial,F);
      FitOpt.IterationPrintFunction(str);
    end
  end
    
  info.bestx = bestx;
  info.minF = Fmin;
  info.nEvals = iTrial;
  info.iter = iTrial;
  info.newbest = newbest;
  if ~isempty(FitOpt.IterFcn)
    UserStop = FitOpt.IterFcn(info);
  else
    UserStop = false;
  end
  if UserStop, stopCode = 2; end
  
  elapsedTime = (cputime-startTime)/60;
  if elapsedTime>FitOpt.maxTime, stopCode = 1; end
  if F<FitOpt.TolFun, stopCode = 3; end
  
  if stopCode, break; end
  
end

if FitOpt.Verbosity>0
  switch stopCode
    case 0, msg = 'Terminated: maximum number of trials reached.';
    case 1, msg = sprintf('Terminated: Time limit of %f minutes reached.',FitOpt.maxTime);
    case 2, msg = 'Terminated: Stopped by user.';
    case 3, msg = sprintf('Terminated: Found a parameter set with error less than %g.',FitOpt.TolFun);
  end
  FitOpt.InfoPrintFunction(msg);
end

info.nTrials = iTrial;
info.elapsedTime = elapsedTime;

end
