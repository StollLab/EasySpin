function ok = test()

% Full D tensor gives the same result as principal values plus DFrame

Sys.S = 1;
Sys.g = [2 2.01 2.02];
Sys.D = [300 50];  % MHz
Sys.DFrame = [0.2 0.5 0.9];

Exp.Field = 350;
Exp.SampleFrame = [0 0 0; 0.3 0.7 1.1];

Pos1 = resfreqs_perturb(Sys,Exp);

R_D2M = erot(Sys.DFrame).';
Dpv = Sys.D(1)*[-1/3 -1/3 2/3] + Sys.D(2)*[1 -1 0];
Sys.D = R_D2M*diag(Dpv)*R_D2M.';
Sys.DFrame = [0 0 0];
Pos2 = resfreqs_perturb(Sys,Exp);

ok = areequal(Pos1,Pos2,1e-8,'rel');
