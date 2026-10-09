function ok = test()

% Full D tensor (also with non-zero trace) is supported by perturbation theory

clear
Sys.S = 1;
Sys.Nucs = '1H';
Sys.g = 2;
Sys.A = [5 6 20];  % MHz
Sys.D = [300 40];  % MHz
Sys.DFrame = [0.2 0.5 0.9];
Exp.Field = 3400;
Exp.mwFreq = 95;
Exp.SampleFrame = [0.3 0.7 1.1];

% principal values + frame
a = endorfrq_perturb(Sys,Exp);

% full D tensor
R_D2M = erot(Sys.DFrame).';
Dpv = Sys.D(1)*[-1/3 -1/3 2/3] + Sys.D(2)*[1 -1 0];
Sys.D = R_D2M*diag(Dpv)*R_D2M.';
Sys.DFrame = [0 0 0];
b = endorfrq_perturb(Sys,Exp);
ok(1) = areequal(sort(a),sort(b),1e-8,'rel');

% full D tensor with non-zero trace: isotropic part has no effect
Sys.D = Sys.D + 100*eye(3);
c = endorfrq_perturb(Sys,Exp);
ok(2) = areequal(sort(a),sort(c),1e-8,'rel');
