function [err,data] = test(opt,olddata)

%=======================================================
% Check a case that gave NaN values in Liouville matrix
% due to overflow in 3j calculation in chili_lm0/lm1/lm2
%=======================================================
Sys.g = [2.008,2.003];
Sys.lw = 0.01;
Sys.tcorr = 10e-6;

Exp.mwFreq = 9.5;
Exp.CenterSweep = [338.4 1.5];

Opt.LLMK = [50 1 1 1];

[x1,y1] = chili(Sys,Exp,Opt);
Sys.A = 10; Sys.Nucs = '1H';
[x1,y1] = chili(Sys,Exp,Opt);
Sys.A = [30 12]; Sys.Nucs = '1H,1H';
[x1,y1] = chili(Sys,Exp,Opt);

err = 0;
data = [];
