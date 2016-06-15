% Basic example of spectral fitting using EasySpin
%====================================================================
clear

% Since we don't have any experimental data available,
% let's create some fake data by simulating a spectrum
% and adding some noise. If you use this example as a
% basis for your fittings, replace this code by code
% that loads your experimental data.

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [120 50 78];
Sys.lwpp = 1;
Exp.mwFreq = 10;
Exp.Range = [300 380];
[B,spc] = pepper(Sys,Exp);
spc = addnoise(spc,50,'n');

% Now we set up the least-squares fitting.
% First comes a starting set of parameters (which we
% obtain by copying the spin system from the simulation
% and changing a few values)
Sys0 = Sys;
Sys0.A = [100 50 78];
Sys0.g(1) = 1.98;

% Next, we specify which parameter we want to be fitted
% and by how much the fitting algorithm can vary it approximately.
Vary.g = [0.1 0 0];
Vary.A = [30 0 0];
SimOpt.Method = 'perturb';

% Calling the fitting function
FitOpt.Method = 'simplex int'; % simplex algorithm, integrals of spectra
esfit('pepper',spc,Sys0,Vary,Exp,SimOpt,FitOpt);
