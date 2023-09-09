function ok = test()

% Test whether pepper returns the correct number of subspectra
% when Opt.separate='orientations'. This only works for crystals.

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
Exp.MolFrame = [20 70 130]*pi/180;
Exp.CrystalSymmetry = 100;

% Two crystal orientations
Exp.SampleFrame = [0 0 0; 234 131 59]*pi/180;

Opt.separate = 'orientations';
Opt.Method = 'perturb2';

[B,spc] = pepper(Sys1,Exp,Opt);
ok(1) = size(spc,1)==2;  % 1 component, 2 orientations each

[B,spc] = pepper({Sys1,Sys3},Exp,Opt);
ok(2) = size(spc,1)==4;  % 2 components, 2 orientations each

[B,spc] = pepper(Sys2,Exp,Opt);
ok(3) = size(spc,1)==8;  % 4 components, 2 orientations each

[B,spc] = pepper({Sys1,Sys2},Exp,Opt);
ok(4) = size(spc,1)==10;  % 4+1 components, 2 orientations each

end
