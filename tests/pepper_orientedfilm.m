function ok = test()

% Test whether sample rotation works correctly for a sample with
% partial ordering.

Sys.S = 1/2;
Sys.lw = 1;

Exp.Ordering = 2;
Exp.mwFreq = 9.5;
Exp.Harmonic = 0;

Opt.GridSize = 31;

Sys.g = [2 2.1 2.2];
Exp.SampleRotation = 0;
[B,spc1] = pepper(Sys,Exp,Opt);

Sys.g = Sys.g([1 3 2]);
Exp.SampleRotation = pi/2;
[B,spc2] = pepper(Sys,Exp);

ok = areequal(spc1,spc2,1e-2,'rel');
