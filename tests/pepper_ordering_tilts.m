function ok = test()

% Test whether sample rotation for a sample with partial ordering works
% consistently with Exp.SampleFrame and Exp.SampleRotation.

Sys.S = 1/2;
Sys.lw = 1;

Exp.Ordering = 3;
Exp.mwFreq = 9.5;
Exp.Harmonic = 0;

Opt.GridSize = 11;
Opt.GridSymmetry = 'Ci';

Sys.g = [2 2.1 2.2];
Sys.gFrame = deg2rad([10 40 20]);

rho = deg2rad(37);  % sample rotation angle
n = [1 4 7];  % sample rotation axis

Exp.SampleRotation = {n,rho};
Exp.SampleFrame = [];
[B,spc1] = pepper(Sys,Exp,Opt);

ang = eulang(rotaxi2mat(n,rho));

Exp.SampleRotation = [];
Exp.SampleFrame = ang;
[B,spc2] = pepper(Sys,Exp,Opt);

ok = areequal(spc1,spc2,1e-2,'rel');
