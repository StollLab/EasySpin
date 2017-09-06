clear, clc

% Define a variable to store error in fitting 
global smallestError
smallestError = Inf;

% Number of data points
N = 500;

% Define a function to fit
xdata = linspace(-5, 5, N);

% Simulate noisy data
x = [1; -1];  % need a vector of parameters
ydata = quadratic(xdata, x);
SNR = 10;
ydata = addnoise(ydata,SNR,'n');

% Define variable ranges for fitting
Vary{1} = [0, 2];  % used for x(1)
Vary{2} = [-4, 0];  % used for x(2)

% Initial guess
x0 = [2; -3];

% esfit algorithms accept and return values that are scaled to the domain 
% [-1,1], so transform them to the [-1,1] domain before and back to their
% original domains afterward
x0 = transformVars(x0, Vary, 0);

% Set some fitting options
FitOpt = struct('nTrials', 10000, 'TolFun', 1e-1);

% Regarding the function to be evaluated at each step in fitting (errorfun 
% in this case), for any input arguments other than the fitting parameters 
% (xdata, ydata, Vary in this case), each esfit algorithm also requires
% them to be input arguments, which can be stored in a cell for convenience
errorfunArgs = {Vary, ydata, xdata};

% Perform fitting
nParameters = numel(x);  % Needed for esfit_montecarlo
xfit = esfit_montecarlo(@errorfun, nParameters, FitOpt, errorfunArgs{:});

% Transform fit parameters back to their original domains
xfit = transformVars(xfit, Vary, 1);

% Plotting

yfit = quadratic(xdata, xfit);

subplot(2,1,1)
hold all
data = plot(xdata, ydata, 'Color', 'black');
fit = plot(xdata, yfit, 'Color', 'red');
legend([data, fit], 'Data', 'Fit')
xlabel('x')
hold off

subplot(2,1,2)
plot(xdata, ydata-yfit)
ylabel('Residuals')
xlabel('x')


% Helper functions
% -------------------------------------------------------------------------

function y = quadratic(xdata, params)
% function to simulate data and fit

y = params(1) + params(2).*xdata.^2;

end

function varargout = errorfun(x, Vary, ydata, xdata)
% this function is fed to the esfit algorithm for checking the error during 
% the fitting process
% Note: the first three arguments (x, Vary, ydata) are always required by 
% the esfit algorithms, whereas any additional arguments are dependent on 
% whether the user-supplied function requires them

global smallestError

xp = transformVars(x, Vary, 1);

simdata = quadratic(xdata, xp);  % the righthand side of this line is the 
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

function xp = transformVars(x, Vary, mode)
% mode 0: transform variables from Vary domains to [-1,1] domain
% mode 1: transform variables from [-1,1] domain to Vary domains

nVars = numel(Vary);

VaryCents = zeros(nVars, 1);  % centers of Vary domains
VaryRanges = zeros(nVars, 1);  % ranges of Vary domains

for iVar=1:nVars
  VaryCents(iVar) = mean(Vary{iVar});
  VaryRanges(iVar) = findRange(Vary{iVar});
end

if mode
  xp = x.*VaryRanges + VaryCents;
else
  xp = (x - VaryCents)./VaryRanges;
end

end

function range = findRange(x)

range = (max(x) - min(x))/2;

end