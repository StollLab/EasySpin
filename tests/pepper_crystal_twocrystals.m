function ok = test()

% Two crystals with Th symmetry and rhombic g

clear Sys Exp Opt
g = [2.0 2.1 2.2];
gFrame = [30 40 78]*pi/180;

Sys.g = g;
Sys.lw = 0.5;
Exp.mwFreq = 9.5;
Exp.Range = [290 350];

Exp.CrystalOrientation = rand(2,3)*pi;
Exp.CrystalSymmetry = 130;

Sys.gFrame = gFrame;
Exp.MolFrame = [0 0 0];
[x,y] = pepper(Sys,Exp);

ok = true;
