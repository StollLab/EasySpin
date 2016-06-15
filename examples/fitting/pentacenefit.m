% fitting the spectrum of pentacene anion radical

% Here is an example you can play with: the simultaneous fitting
% of four hyperfine constants.
% It represents a rather difficult fitting problem,
% if no pre-analysis of the spectrum is done to narrow
% down the ranges of the hyperfine coupling constants.
% Bottom line: always analyse your experimental spectra
% as thoroughly as possible before you try to least-squares
% fit the spectra.

% Compute a noisy spectrum
clear
Sys = struct('g',2,'Nucs','1H,1H,1H,1H','n',[2 4 4 4]);
Sys.lwpp = [0,0.01]; % Lorentzian lines
Sys.A =  [11.9 8.5 2.6 2.4];
Exp = struct('mwFreq',9.5,'CenterSweep',[339.4 4],'nPoints',2e3);
Exp.Harmonic = 1;
[xa,ya] = garlic(Sys,Exp);
ya = addnoise(ya,150,'n');

% Modify the hyperfine values and set up a rather wide
% search range for them.

Sys.A = [12 8.5 2.6 2.4];
Vary.A = [2 2 1 1];

FitOpt.Method = 'genetic int';
BestSys = esfit('garlic',ya,Sys,Vary,Exp,[],FitOpt);
