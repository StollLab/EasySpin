% Effect of ordering potential on slow-motion an rigid-limit EPR spectra
%===============================================================================
% This example illustrates how the ordering potential input works and how an
% ordering potential affects a simple EPR spectrum, both in the rigid limit and
% in the slow-motion regime.

clear, clc, clf

% Define a simple spin system
Sys.g = [2 2.05];
Sys.lw = [0 0.5];

Exp.mwFreq = 9.5; % GHz
Exp.Range = [320 350]; % mT
Exp.CrystalOrientation = [0 0 0];
Exp.Harmonic = 0;   % show absorption spectrum, which is more intuitive

Opt.LLMK = [100 0 0 0]; % orientational basis set

% Set potential coefficient for L=2, M=2, K=2 term
c200 = 2; % positive indicates preferential along z axis

% Simulate reference rigid-limit spectrum without ordering
Exp.Ordering = 0;
[B,spcref] = pepper(Sys,Exp,Opt);

% Simulate rigid-limit spectrum
Exp.Ordering = c200;
[B,spc0] = pepper(Sys,Exp,Opt);

% Simulate quasi-rigid spectrum using chili
Sys.Potential = [2 0 0 c200];
Sys.logtcorr = -5; % log10(tcorr/seconds)
[B,spc1] = chili(Sys,Exp,Opt);

% Simulate slow-motion spectrum using chili
Sys.Potential = [2 0 0 c200];
Sys.logtcorr = -9;
[B,spc2] = chili(Sys,Exp,Opt);

plot(B,spcref,B,spc0,B,spc1,B,spc2);
legend('no orderiing','rigid limit (pepper)','quasi-rigid (chili)','slow motion (chili');
xlabel('magnetic field (mT)');
