% Fitting the spectrum vs fitting its integral

clear, clc

% This example illustrates how important it is to use good
% starting values and to tell esfit whether to fit the
% spectrum directly or to fit its integral.

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
Sys0.g(3) = 2.21;
Vary.g = [0 0 0.04];
SimOpt.Method = 'perturb';

% Using 'fcn', the spectra are directly compared. For
% 'int', the spectra are integrated before comparison.
% 'fcn' ends up in wrong minima much more often than 'int'.

FitOpt.Method = 'simplex fcn';
esfit(y,@pepper,{Sys0,Exp,SimOpt},{Vary},FitOpt);

FitOpt.Method = 'simplex int';
esfit(y,@pepper,{Sys0,Exp,SimOpt},{Vary},FitOpt);
