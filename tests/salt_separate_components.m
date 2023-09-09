function ok = test()

% Test whether salt returns the correct number of subspectra
% when Opt.separate='components'.

% Define several spin systems
Sys1.Nucs = '35Cl';
Sys1.A = [3 6];
Sys1.lwEndor = 0.2;

Sys2.Nucs = '37Cl';
Sys2.A = 5;
Sys2.lwEndor = 0.2;

Sys3.Nucs = 'Cl,Cl';
Sys3.A = [2 4];
Sys3.lwEndor = 0.2;

% Experimental settings
Exp.Field = 3400;  % mT
Exp.Range = [0 40];  % MHz

Opt.Method = 'perturb2';

Opt.separate = 'components';

% Powder
[nu,spc] = salt(Sys1,Exp,Opt);
ok(1) = size(spc,1)==1;

[nu,spc] = salt({Sys1,Sys2},Exp,Opt);
ok(2) = size(spc,1)==2;

[nu,spc] = salt(Sys3,Exp,Opt);
ok(3) = size(spc,1)==4;

[nu,spc] = salt({Sys1,Sys3},Exp,Opt);
ok(4) = size(spc,1)==5;

% Crystal
Exp.SampleFrame = [20 70 130]*pi/180;
[nu,spc] = salt(Sys1,Exp,Opt);
ok(5) = size(spc,1)==1;

[nu,spc] = salt({Sys1,Sys2},Exp,Opt);
ok(6) = size(spc,1)==2;

[nu,spc] = salt(Sys3,Exp,Opt);
ok(7) = size(spc,1)==4;

[nu,spc] = salt({Sys1,Sys3},Exp,Opt);
ok(8) = size(spc,1)==5;

end
