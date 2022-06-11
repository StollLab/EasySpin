% Basic example of spectral fitting using EasySpin
%====================================================================
clear, clc

% In place of loading experimental data, we create some mock data by
% simulating a spectrum and adding some noise. If you use this example
% as a basis for your fittings, replace this code by code that loads your
% experimental data.

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [120 50 78];  % MHz
Sys.lwpp = 1;  % mT
Exp.mwFreq = 10;  % GHz
Exp.Range = [300 380];  % mT
[B,spc] = pepper(Sys,Exp);
spc = addnoise(spc,150,'n');

% Now we set up the least-squares fitting. First we specify a spin system
% that is going to be used as a starting point of the fit.
Sys0.g = [1.98 2.1 2.2];
Sys0.Nucs = '1H';
Sys0.A = [100 50 78];  % MHz
Sys0.lwpp = 1;  % mT

% Next, we specify a second spin system structure that contains only the
% parameters we want to be fitted and by how much they can be varied
% relative to the value in Sys0.
SysVary.g(1) = 0.1;
SysVary.A = [30 0 0];  % MHz

% We also provide options for the simulation function.
SimOpt.Method = 'perturb';

% Finally, we specify the fitting algorithm and what should be fitted.
FitOpt.Method = 'simplex int'; % simplex algorithm, integral of spectrum

% Calling esfit without outputs will start the GUI.
esfit(spc,@pepper,{Sys0,Exp,SimOpt},{SysVary},FitOpt);

% If the fitting algorithm doesn't find the correct minimum, you can change
% the algorithm, target function, and starting point in the UI. For example,
% run Monte Carlo for a bit, save the best fit, and then use that as the
% starting point with Nelder-Mead to close in on the correct fit.

% As an alternative to the GUI, call esfit with one output. Then esfit will
% run without GUI and return the fit result in the output.
fit = esfit(spc,@pepper,{Sys0,Exp,SimOpt},{SysVary},FitOpt);
