function ok = test()

% Test whether salt returns the correct number of subspectra
% when Opt.Output='orientations'. This only works for crystals.

% Define several spin systems
Sys1.Nucs = '35Cl';
Sys1.A = [3 6];
Sys1.lwEndor = 0.2;

Sys2.Nucs = 'Cl,Cl';
Sys2.A = [2 4];
Sys2.lwEndor = 0.2;

Sys3.Nucs = '37Cl';
Sys3.A = 5;
Sys3.lwEndor = 0.2;

% Experimental settings
Exp.Field = 3400;  % mT
Exp.Range = [0 20];
Exp.MolFrame = [20 70 130]*pi/180;
Exp.CrystalSymmetry = 100;

% Two crystal orientations
Exp.SampleFrame = [0 0 0; 234 131 59]*pi/180;

Opt.Output = 'orientations';
Opt.Method = 'perturb2';

[nu,spc] = salt(Sys1,Exp,Opt);
ok(1) = size(spc,1)==2;  % 1 component, 2 orientations each

[nu,spc] = salt({Sys1,Sys3},Exp,Opt);
ok(2) = size(spc,1)==4;  % 2 components, 2 orientations each

[nu,spc] = salt(Sys2,Exp,Opt);
ok(3) = size(spc,1)==8;  % 4 components, 2 orientations each

[nu,spc] = salt({Sys1,Sys2},Exp,Opt);
ok(4) = size(spc,1)==10;  % 4+1 components, 2 orientations each

end
