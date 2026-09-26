function ok = test()

% Thermal polarization for S>1/2 agrees with resfreqs_matrix

Sys.S = 3/2;
Sys.g = 2;
Sys.D = 30;  % MHz

Exp.Field = 350;  % mT
Exp.Range = [9 10.5];  % GHz
Exp.SampleFrame = [0 0 0];
Exp.Temperature = 2;  % K

[Pos_p,Int_p] = resfreqs_perturb(Sys,Exp);
[Pos_m,Int_m] = resfreqs_matrix(Sys,Exp);

[Pos_p,idx] = sort(Pos_p(:)); Int_p = Int_p(idx);
[Pos_m,idx] = sort(Pos_m(:)); Int_m = Int_m(idx);

ok(1) = areequal(Pos_p,Pos_m,1e-3,'rel');
ok(2) = areequal(Int_p/max(Int_p),Int_m/max(Int_m),2e-2,'rel');  % perturbation theory neglects D in populations and transition rates
