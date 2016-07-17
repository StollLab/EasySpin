% Fitting the spectrum vs fitting its integral

clear, clf

% This example illustrates how important it is to tell esfit
% whether to fit the spectrum directly or to fit its
% integral.

% Generate a simple spectrum of one unpaired electron and
% a hydrogen nucleus
Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [50 100 150];
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.CenterSweep = [330 150];
[x,y] = pepper(Sys,Exp);
y = addnoise(y,60,'n');

% For fitting, slightly shift the g value.

Sys0 = Sys;
Sys0.g(3) = 2.22;
Vary.g = [0 0 0.02];
SimOpt.Method = 'perturb';

% Using 'fcn', the spectra are directly compared. For
% 'int', the spectra are integrated before comparison.
% 'fcn' ends up in wrong minima much more often than 'int'.

FitOpt.Method = 'simplex fcn';
esfit('pepper',y,Sys0,Vary,Exp,SimOpt,FitOpt);

FitOpt.Method = 'simplex int';
esfit('pepper',y,Sys0,Vary,Exp,SimOpt,FitOpt);
