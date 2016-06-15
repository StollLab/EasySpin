function [err,data] = test(opt,olddata)

% Test whether Opt.Output = 'separate' works with all spectrum
% construction methods.

Sys.Nucs = '1H';
Sys.A = 10;
Sys.lwpp = [0 0.01];
Exp.mwFreq = 9.5;
Exp.Range = [337.5 340.5];

Opt.Output = 'separate';

Opt.AccumMethod = 'explicit';
[x,y0] = garlic(Sys,Exp,Opt);
Opt.AccumMethod = 'template';
[x,y1] = garlic(Sys,Exp,Opt);
Opt.AccumMethod = 'binning';
[x,y2] = garlic(Sys,Exp,Opt);

y0s = sum(y0);
y1s = sum(y1);
y2s = sum(y2);

err = false;
data = [];

