% Using esfit for non-EasySpin model functions
%===============================================================================
% This example demonstrates how to use esfit to fit an arbitrary
% (non-EasySpin) model function that takes a single parameter vector
% as input.

clear, clc, clf

% Define a model function
%-------------------------------------------------------------------------------
t = linspace(0, 6, 500);
model = @(p) exp(-(t/p(1)).^p(2));

% Generate synthetic data
%-------------------------------------------------------------------------------
params = [3, 2.7];   % parameters
ydata = model(params);
ydata = addnoise(ydata,60,'n');  % add some noise

% Set parameter
%-------------------------------------------------------------------------------
% Initial parameter guess
params0 = [1, 2];

% Define parameter ranges (lower and upper bounds) for fitting
params_lb = [0.1 1];
params_ub = [6 3];

% Perform fitting
result = esfit(ydata, model, params0, params_lb, params_ub);

% Plotting
subplot(3,1,[1 2])
plot(t, ydata, t, result.fit);
legend('data', 'fit')
legend boxoff

subplot(3,1,3)
plot(t, ydata(:)-result.fit(:))
ylabel('Residuals')
