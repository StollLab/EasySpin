function ok = test()

% Test whether pepper returns the correct number of subspectra
% when Opt.Output='transitions'. This only works for powders.

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

Opt.Output = 'transitions';
Opt.Method = 'perturb2';

%35Cl and 37Cl both have spin 3/2, so 4 allowed EPR transitions

% Powder
[B,spc] = pepper(Sys1,Exp,Opt);
ok(1) = size(spc,1)==4;  % 1 components,

[B,spc] = pepper({Sys1,Sys3},Exp,Opt);
ok(2) = size(spc,1)==5;

[B,spc] = pepper(Sys2,Exp,Opt);
ok(3) = size(spc,1)==64;  % 4 isotopologue components, 4x4=16 transitions each

[B,spc] = pepper({Sys1,Sys2},Exp,Opt);
ok(4) = size(spc,1)==68;  % 4 components with 16 transitions, plus 1 comp. with 4

end
