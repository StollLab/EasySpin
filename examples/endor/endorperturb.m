% ENDOR simulations: matrix diagonalization vs perturbation theory
%======================================================================
% For ENDOR spectra with many nuclei, e.g. protons, ENDOR simulations
% can become very slow. In this case, you can try to use perturbation
% theory instead of the default matrix diagonalization to speed up the
% simulation.
% Perturbation theory yields correct results as long as the hyperfine
% coupling is much smaller than the microwave frequency.

clear, clf, clc

% Spin system parameters
Sys.g = [2.3,2.1,2];
Sys = nucspinadd(Sys,'1H',[3,6,2]);
Sys = nucspinadd(Sys,'1H',[5,3,1]);
Sys = nucspinadd(Sys,'1H',[8,6,-4]);
Sys = nucspinadd(Sys,'1H',[1,-1,1.6]);
Sys.HStrain = [1 1 1]*200;     % MHz
Sys.lwEndor = 0.1;    % MHz

% ENDOR parameters
Exp.Field = 300;
Exp.Range = larmorfrq('1H',Exp.Field) + [-1 1]*10; % MHz
Exp.mwFreq = 9.5;
Exp.ExciteWidth = 100;

% Simulation options
Opt.nKnots = [31 0]; % higher orientational resolution

% (1) matrix diagonalization: slow
Opt.Method = 'matrix';
tic
[freq,spec_m] = salt(Sys,Exp,Opt);
toc

% (2) perturbation theory: faster
Opt.Method = 'perturb1';
tic
[freq,spec_p] = salt(Sys,Exp,Opt);
toc

% Plotting
freq = freq - larmorfrq('1H',Exp.Field);
plot(freq,spec_m/max(spec_m),'b',freq,spec_p/max(spec_p),'g');

xlabel('frequency [MHz]');
title('ENDOR spectra');
legend('matrix diagonalization','perturbation theory');
