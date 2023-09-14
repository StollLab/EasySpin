% Fitting magnetic susceptibility curve
%===============================================================================

clear, clc

% Generate some data
dataSys.S = 3/2;
dataSys.D = -20*clight*1e-4;
dataSys.g = 2.3;

% Temperatures and field values
T = [2:0.5:20, 25:5:300];
Exp.Temperature = T;
Exp.Field = [100, 5000];

% Calculate molar susceptibility
Opt.Output = 'chimolT';
Opt.Units = 'CGS';
chimolT  = curry(dataSys,Exp,Opt);

% Add some noise
chimolT = addnoise(chimolT,80,'n');

% Plot the data
plot(T,chimolT,'.')
xlabel('temperature (K)')
ylabel('\chi T  (K cm^3 mol^{-1})')
grid on

%% 
% Make an initial guess for the parameters
Sys0 = dataSys;
Sys0.S = 3/2;
Sys0.D = -30*clight*1e-4;
Sys0.g = 2;

% What should be varied during the fit
vSys.D = 40*clight*1e-4;
vSys.g = 0.5;

% No autoscaling! the absolute value of the data contains information
FitOpt.AutoScale = 'none';

% Run least-squares fitting
fitresult = esfit(chimolT(:),@curry,{Sys0,Exp,Opt},{vSys},FitOpt)

%%
% Plotting
fit = reshape(fitresult.fit,2,[]);
plot(T,chimolT,'.',T,fit)
xlabel('temperature (K)')
ylabel('\chi T  (K cm^3 mol^{-1})')
grid on
