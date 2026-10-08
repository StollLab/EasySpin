function ok = test()

% H strain

Sys.g = [2 2.1 2.2];
Sys.HStrain = [10 20 30];  % MHz

Exp.Field = 350;
Exp.SampleFrame = [0 0 0];  % B0 along molecular z

[~,~,Wid] = resfreqs_perturb(Sys,Exp);
ok(1) = areequal(Wid,30,1e-10,'rel');

