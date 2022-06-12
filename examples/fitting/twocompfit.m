% Fitting a two-component spectrum
%====================================================================
clear, clc

% Since we don't have any experimental data available, we create a
% two-component spectrum and add noise. If you use this example as
% a basis for your fittings, replace this code by code that loads
% your experimental data.

Sys1.g = [2 2.1 2.2];
Sys1.lwpp = 1;  % mT
Sys1.weight = 0.7;
Sys2.g = [1.85 1.9 2.05];
Sys2.lwpp = 2;  % mT
Sys2.weight = 0.3;

Exp.mwFreq = 10;  % GHz
Exp.Range = [300 400];  % mT

[B,spc] = pepper({Sys1,Sys2},Exp);
spc = addnoise(spc,50,'n');
plot(B,spc);

% Next we set up the least-squares fitting.
% First comes a starting set of parameters (which we obtain by reusing
% the spin systems from the simulation and changing a few values)
Sys1.g(1) = 1.98;
Sys2.lwpp = 1.3;
Sys2.weight = 0.6;
Sys = {Sys1,Sys2};

% Next, we specify which parameter we want to be fitted and by how much
% the fitting algorithm can vary it.
Vary1.g = [0.03 0 0];
Vary2.lwpp = 0.9;
Vary2.weight = 0.3;
Vary = {Vary1,Vary2};

% Call the fitting function
SimOpt.Method = 'perturb';
FitOpt.Method = 'simplex int'; % simplex algorithm, integrals of spectra
result = esfit(spc,@pepper,{Sys,Exp,SimOpt},{Vary},FitOpt);

% Plot the fit result
plot(B,spc,B,result.fit)
