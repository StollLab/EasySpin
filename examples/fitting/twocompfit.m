% Fitting a two-component spectrum
%====================================================================
clear

% Since we don't have any experimental data available, we create a
% two-component spectrum and add noise. If you use this example as
% a basis for your fittings, replace this code by code that loads
% your experimental data.

Sys1.g = [2 2.1 2.2];
Sys1.lwpp = 1;
Sys1.weight = 0.7;
Sys2.g = [1.85 1.9 2.05];
Sys2.lwpp = 2;
Sys1.weight = 0.3;
Exp.mwFreq = 10;
Exp.Range = [300 400];
[B,spc] = pepper({Sys1,Sys2},Exp);
spc = addnoise(spc,50,'n');
plot(B,spc);

% Next we set up the least-squares fitting.
% First comes a starting set of parameters (which we obtain by reusing
% the spin systems from the simulation and changing a few values)
Sys1.g(1) = 1.98;
Sys2.lwpp = 1;
Sys2.weight = 0.6;

% Next, we specify which parameter we want to be fitted and by how much
% the fitting algorithm can vary it.
Vary1.g = [0.03 0 0];
Vary2.lwpp = 1;
Vary2.weight = 0.3;

% Calling the fitting function
SimOpt.Method = 'perturb';
FitOpt.Method = 'simplex int'; % simplex algorithm, integrals of spectra
esfit('pepper',spc,{Sys1,Sys2},{Vary1,Vary2},Exp,SimOpt,FitOpt);
