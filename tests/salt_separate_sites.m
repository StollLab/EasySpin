function ok = test()

% Test whether salt returns the correct number of subspectra
% when Opt.separate='sites'. This works only for crystals.

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
Exp.CrystalSymmetry = 100; % Space group 100 has 8 sites.

% Two crystal orientations
Exp.SampleFrame = [0 0 0; 234 131 59]*pi/180;

Opt.separate = 'sites';

[B,spc] = salt(Sys1,Exp,Opt);
ok(1) = size(spc,1)==8;  % 1 component with 8 sites

[B,spc] = salt({Sys1,Sys3},Exp,Opt);
ok(2) = size(spc,1)==16;  % 2 components with 8 sites

[B,spc] = salt(Sys2,Exp,Opt);
ok(3) = size(spc,1)==32;  % 2x2=4 components with 8 sites

[B,spc] = salt({Sys1,Sys2},Exp,Opt);
ok(4) = size(spc,1)==40;  % 1+2x2=5 components with 8 sites

end
