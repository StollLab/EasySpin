function ok = test()

% Two crystals with Th symmetry and rhombic g

g = [2.0 2.1 2.2];
gFrame = [30 40 78]*pi/180;

Sys.g = g;
Sys.gFrame = gFrame;
Sys.lw = 0.5;

Exp.SampleFrame = [10 24 61; 222 55 99]*pi/180;
Exp.CrystalSymmetry = 130;

Exp.mwFreq = 9.5;
Exp.Range = [290 350];

[B,spc] = pepper(Sys,Exp);
ok(1) = size(spc,1)==1;

Opt.separate = 'orientations';
[B,spc] = pepper(Sys,Exp,Opt);
ok(2) = size(spc,1)==2;
