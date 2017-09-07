% Utilizing low-level fitting functions directly for custom fitting
%===============================================================================
% This example demonstrates how to use the fitting functions underlying
% EasySpin's esfit (such as esfit_montecarlo etc) to perform least-squares
% fitting.

clear, clc, clf

% Define a variable to store error in fitting 
global smallestError
smallestError = Inf;

% Define a function to fit
%-------------------------------------------------------------------------------
xdata = linspace(-5, 5, 500);  % x axis
params = [1; -1];   % parameters
ydata = myfunction(xdata, params); % custom function, defined below
ydata = addnoise(ydata,10,'n');  % add some noise


% Set up fitting
%-------------------------------------------------------------------------------
% Initial parameter guess
params0 = [1.5; -3];

% Define parameter ranges (min, max) for fitting
Vary{1} = [0, 2];  % for params(1)
Vary{2} = [-4, 0];  % for params(2)

% esfit algorithms accept and return values that are scaled to the domain 
% [-1,1], so transform them to the [-1,1] domain before and back to their
% original domains afterward. The function transformVars is defined below.
params0 = transformVars(params0, Vary, 0);

% Set some fitting options
FitOpt = struct('nTrials', 10000, 'TolFun', 1e-4, 'PrintLevel', 2);

% Regarding the function to be evaluated at each step in fitting (errorfun 
% in this case), for any input arguments other than the fitting parameters 
% (xdata, ydata, Vary in this case), each esfit algorithm also requires
% them to be input arguments, which can be stored in a cell for convenience.
errorfunArgs = {Vary, ydata, xdata};

% Perform fitting
%nParameters = numel(params);  % Needed for esfit_montecarlo, others need params0
paramsfit = esfit_simplex(@errorfun, params0, FitOpt, errorfunArgs{:});

% Transform fit parameters back to their original domains
paramsfit = transformVars(paramsfit, Vary, 1);

% Plotting
%-------------------------------------------------------------------------------
yfit = myfunction(xdata, paramsfit);

subplot(3,1,[1 2])
plot(xdata, ydata, xdata, yfit);
legend('data', 'fit')
legend boxoff
xlabel('x')
hold off

subplot(3,1,3)
plot(xdata, ydata-yfit)
ylabel('Residuals')
xlabel('x')


%===============================================================================
% Helper functions
%===============================================================================

function y = myfunction(xdata, params)
% function to simulate data and fit

y = params(1) + params(2).*xdata.^2;

end

%-------------------------------------------------------------------------------
function varargout = errorfun(x, Vary, ydata, xdata)
% This function is fed to the esfit algorithm for checking the error during 
% the fitting process.
% Note: the first three arguments (x, Vary, ydata) are always required by 
% the esfit algorithms, whereas any additional arguments are dependent on 
% whether the user-supplied function requires them

global smallestError

xp = transformVars(x, Vary, 1);

simdata = myfunction(xdata, xp);  % the righthand side of this line is the 
                                 % only line that can be changed, which 
                                 % depends on what function is being fitted
                                 % and what arguments it requires
resid = ydata - simdata;
rmsd = real(sqrt(mean(resid.^2)));

isNewBest = rmsd < smallestError;

if isNewBest
  smallestError = rmsd;
end

out = {rmsd, resid, simdata};

% some esfit algorithms don't need all of these output arguments
varargout = out(1:nargout);

end

%-------------------------------------------------------------------------------
function xp = transformVars(x, Vary, mode)
% mode 0: transform variables from Vary domains to [-1,1] domain
% mode 1: transform variables from [-1,1] domain to Vary domains

nVars = numel(Vary);

VaryCents = zeros(nVars, 1);  % centers of Vary domains
VaryRanges = zeros(nVars, 1);  % ranges of Vary domains

for iVar = 1:nVars
  v = Vary{iVar};
  VaryCents(iVar) = mean(v);
  VaryRanges(iVar) = (max(v)-min(v))/2;
end

if mode
  xp = x.*VaryRanges + VaryCents;
else
  xp = (x - VaryCents)./VaryRanges;
end

end
