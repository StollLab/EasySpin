function ok = test()

nOrientations = 7;

Sys.g = [2 2.1 2.2];
Sys.lwpp = 0.1;

Exp.MolFrame = [8 20 76]*pi/180;
Exp.CrystalSymmetry = 1;

Exp.SampleFrame = [10 33 -8]*pi/180;
rho = linspace(0,pi,nOrientations);
Exp.SampleRotation = {'x',rho};

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.nPoints = 1e4;

Opt.separate = 'orientations';
[~,spc_separate] = pepper(Sys,Exp,Opt);
Opt.separate = '';
[~,spc] = pepper(Sys,Exp,Opt);

ok(1) = size(spc,1)==1;
ok(2) = size(spc_separate,1)==nOrientations;
