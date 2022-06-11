% fitting the spectrum of pentacene anion radical

% Here is an example you can play with: the simultaneous fitting
% of four hyperfine constants.
% It represents a rather difficult fitting problem,
% if no pre-analysis of the spectrum is done to narrow
% down the ranges of the hyperfine coupling constants.
% Bottom line: always analyse your experimental spectra
% as thoroughly as possible before you try to least-squares
% fit the spectra.

clear, clc

% Compute a noisy spectrum
Pentacene.g = 2;
Pentacene.Nucs = '1H,1H,1H,1H';
Pentacene.n = [2 4 4 4];
Pentacene.lwpp = [0 0.01]; % Lorentzian broadening
Pentacene.A =  [11.9 8.5 2.6 2.4];  % MHz

Exp.mwFreq = 9.5;  % GHz
Exp.CenterSweep = [339.4 5];  % mT
Exp.nPoints = 3e3;
Exp.Harmonic = 1;

[xa,ya] = garlic(Pentacene,Exp);
ya = addnoise(ya,150,'n');

% Modify the hyperfine values and set up a rather wide
% search range for them.

Pentacene.A = [12 8.5 2.6 2.4];
Vary.A = [2 2 1 1];  % MHz

% Run fitting
FitOpt.Method = 'genetic int';
fit = esfit(ya,@garlic,{Pentacene,Exp},{Vary},FitOpt);
