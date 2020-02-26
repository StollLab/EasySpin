function ok = test()

% Test memory reallocation in chili_lm

Sys.g = [2.01 2.02 2.03];
Sys.tcorr = 1e-8;
Sys.lw = 0.1;
Exp.mwFreq = 9.7;

% no nuclei
Opt.AllocationBlockSize = 2e3;
[B,spc] = chili(Sys,Exp,Opt);

% one nucleus
Sys.Nucs = '14N';
Sys.A = [10 20 30];
Opt.AllocationBlockSize = 5e4;
[B,spc] = chili(Sys,Exp,Opt);

% two nuclei
Sys.Nucs = '15N,15N';
Sys.A = [10 10 20; 10 10 20];
Opt.AllocationBlockSize = 5e4;
[B,spc] = chili(Sys,Exp,Opt);

ok = true;
