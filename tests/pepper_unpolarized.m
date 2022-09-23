function ok = test()

% Compare matrix-diag and perturbation-theory simulations
% for excitation with unpolarized light.

Sys.g = [3 2];
Sys.lwpp = 10;

Exp.mwFreq = 9.5;
Exp.Range = [180 380];
Exp.mwMode = {rand*pi/2 'unpolarized'};

Opt.Method = 'matrix';
[x,y1] = pepper(Sys,Exp,Opt);
Opt.Method = 'perturb';
[x,y2] = pepper(Sys,Exp,Opt);

maxy1 = max(y1);
y1 = y1/maxy1;
y2 = y2/maxy1;

ok = max(abs(y2-y1))<1e-6;
