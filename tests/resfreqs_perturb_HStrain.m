function ok = test()

% H strain alone and combined with g strain

Sys.g = [2 2.1 2.2];
Sys.HStrain = [10 20 30];  % MHz

Exp.Field = 350;
Exp.SampleFrame = [0 0 0];  % B0 along molecular z

[~,~,Wid] = resfreqs_perturb(Sys,Exp);
ok(1) = areequal(Wid,30,1e-10,'rel');

Sys.gStrain = [0 0 0.01];
[Pos,~,Wid] = resfreqs_perturb(Sys,Exp);
Wg = Pos*0.01/2.2;  % MHz
ok(2) = areequal(Wid,sqrt(30^2+Wg^2),1e-10,'rel');
