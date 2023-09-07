function ok = test()

% Test whether pepper returns the correct number of subspectra
% when Opt.Output='sites'. This works only for crystals.

% Define several spin systems
Sys1.g = 2;
Sys1.Nucs = '35Cl';
Sys1.A = 10;
Sys1.lwpp = 0.2;

Sys2.g = 2;
Sys2.Nucs = 'Cl,Cl';
Sys2.A = [20 7];
Sys2.lwpp = 0.2;

Sys3.g = 2.02;
Sys3.lwpp = 0.4;

Exp.mwFreq = 9.5;
Exp.Range = [330 360];

Opt.Output = 'sites';

Exp.MolFrame = [20 70 130]*pi/180;
Exp.CrystalSymmetry = 100;

% Space group 100 has 8 sites.

[B,spc] = pepper(Sys1,Exp,Opt);
ok(1) = size(spc,1)==8;  % 1 component with 8 sites

[B,spc] = pepper({Sys1,Sys3},Exp,Opt);
ok(2) = size(spc,1)==16;  % 2 components with 8 sites

[B,spc] = pepper(Sys2,Exp,Opt);
ok(3) = size(spc,1)==32;  % 2x2=4 components with 8 sites

[B,spc] = pepper({Sys1,Sys2},Exp,Opt);
ok(4) = size(spc,1)==40;  % 1+2x2=5 components with 8 sites

end
