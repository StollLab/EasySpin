function ok = test()

% g strain widths along principal axes agree with resfreqs_matrix

Sys.g = [2 2.1 2.2];
Sys.gStrain = [0.01 0.02 0.03];

Exp.Field = 350;
Exp.SampleFrame = [0 pi/2 0; pi/2 pi/2 0; 0 0 0];  % B0 along x, y, z

[Pos_p,~,Wid_p] = resfreqs_perturb(Sys,Exp);
[Pos_m,~,Wid_m] = resfreqs_matrix(Sys,Exp);

ok(1) = areequal(Pos_p,Pos_m,1e-6,'rel');
ok(2) = areequal(Wid_p,Wid_m,1e-6,'rel');
