% Effect of orientational potential on slow-motion an rigid-limit EPR spectra
%===============================================================================
% This example illustrates how the orientational potential input works and how an
% orientational potential affects a simple EPR spectrum, both in the rigid limit and
% in the slow-motion regime.

clear, clc, clf

% Define a simple spin system
Sys.g = [2 2.05];
Sys.lw = [0.5 0];

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [325 345];  % mT
Exp.Harmonic = 0;   % show absorption spectrum, which is more intuitive

% orientational basis set of chili simulations
% (three zeros here due to the specific SampleFrame and the axial g;
% in general non-zero values will be needed!)
Opt.LLMK = [80 0 0 0]; 

% Set potential coefficient for L=2, M=2, K=2 term
c200 = 2; % positive indicates preferential along z axis

% Simulate reference rigid-limit spectrum without orientational potential
Exp.Ordering = 0;
[B,spcref] = pepper(Sys,Exp,Opt);

% Simulate rigid-limit spectrum
Exp.Ordering = c200;
[B,spc0] = pepper(Sys,Exp,Opt);

% Simulate quasi-rigid spectrum using chili
Exp.SampleFrame = [0 0 0];
Sys.Potential = [2 0 0 c200];
Sys.tcorr = 1e-4;  % s
[B,spc1] = chili(Sys,Exp,Opt);

% Simulate slow-motion spectrum using chili
Sys.Potential = [2 0 0 c200];
Sys.tcorr = 3e-9;  % s
[B,spc2] = chili(Sys,Exp,Opt);

plot(B,spcref,B,spc0,B,spc1,B,spc2);
legend('no ordering','ordering; rigid limit (pepper)','ordering; quasi-rigid (chili)','ordering; slow motion (chili)');
xlabel('magnetic field (mT)');
