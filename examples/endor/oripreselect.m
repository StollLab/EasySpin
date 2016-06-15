% Speed-up of ENDOR simulations using orientation pre-selection
%======================================================================
% Opt.OriPreSelect provides a speed-up for simulations of
% strongly orientation-selective ENDOR spectra. If the orientation
% selection is not strong, the speed-up is only minor.

clear, clf, clc

% Spin system
Sys.g = [2.3,2.1,2];
Sys = nucspinadd(Sys,'1H',[3,6,2]);
Sys = nucspinadd(Sys,'1H',[5,3,1]);
Sys = nucspinadd(Sys,'1H',[8,6,-4]);
Sys.HStrain = [1 1 1]*100;     % MHz
Sys.lwEndor = 0.1;    % MHz

% ENDOR simulations
Exp.Field = 297;
Exp.CenterSweep = [larmorfrq('1H',Exp.Field), 12]; % MHz
Exp.mwFreq = 9.5;
Exp.ExciteWidth = 50;

% Simulation options
Opt.nKnots = [61 0];

% (1) no orientation pre-selection: slow
Opt.OriPreSelect = 0;
tic
[x,y1] = salt(Sys,Exp,Opt);
toc

% (2) with orientation pre-selection: faster
Opt.OriPreSelect = 1;
tic
[x,y2] = salt(Sys,Exp,Opt);
toc

% Plotting
plot(x,y1,'g',x,y2,'r')

xlabel('frequency [MHz]');
title('ENDOR spectra');
legend('no pre-selection','pre-selection');
