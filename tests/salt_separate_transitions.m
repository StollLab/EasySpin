function ok = test()

% Test whether salt returns the correct number of subspectra
% when Opt.separate='transitions'. This only works for powders.

% Define several spin systems
Sys1.Nucs = '15N';
Sys1.A = [3 6];
Sys1.lwEndor = 0.2;

Sys2.Nucs = '37Cl';
Sys2.A = [5 7];
Sys2.lwEndor = 0.2;

Sys3.Nucs = 'Cl';
Sys3.A = [2 6];
Sys3.lwEndor = 0.2;

% Experimental settings
Exp.Field = 3400;  % mT
Exp.Range = [0 40];

Opt.separate = 'transitions';
Opt.Method = 'perturb2';

%35Cl and 37Cl both have spin 3/2, so 4 allowed EPR transitions

% Powder
[nu,spc] = salt(Sys1,Exp,Opt);
ok(1) = size(spc,1)==2;  % 1 component, 2 transitions (I=1/2)

[nu,spc] = salt(Sys2,Exp,Opt);
ok(1) = size(spc,1)==6;  % 1 component, 6 transitions (I=3/2)

[nu,spc] = salt({Sys1,Sys2},Exp,Opt);
ok(2) = size(spc,1)==8;  % 2 components, 2+6 transitions

[nu,spc] = salt(Sys3,Exp,Opt);
ok(3) = size(spc,1)==12;  % 2 isotopologue components, 6+6 transitions

[nu,spc] = salt({Sys1,Sys3},Exp,Opt);
ok(4) = size(spc,1)==14;  % 3 components, 6+6+2 transitions

end
