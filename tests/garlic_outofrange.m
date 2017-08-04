function [err,data] = test(opt,olddata)

% Make sure that garlic doesn't crash if there are out-of-range transitions.
% This was a bug in 5.0.16.

Sys.Nucs = '14N,14N,1H,1H,1H';
Sys.n = [2 2 2 2 6];
Sys.g = [2.004];

A_N1 = 3.18; % in MHz
A_N2 = 10.97;
A_H1 = 7.60;
A_H2 = 9.06;
A_H3 = 12.14;

Sys.A = [A_N1;A_N2;A_H1;A_H2;A_H3];

Sys.lwpp = [0.0 0.00156];

Exp.Range = [347 353];
Exp.mwFreq = 9.85184;
Exp.ModAmp = 0.03;

Opt.Method = 'perturb';

[x,y] = garlic(Sys,Exp,Opt);

if opt.Display
  plot(x,y);
end

err = false;
data = [];
