function ok = test()

% Non-symmetric full g and A tensors: perturbation theory agrees with matrix diagonalization

clear
Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.g = [2 0.1 0.05; -0.08 2.1 0.03; 0.02 -0.06 2.2];
Sys.A = [10 30 -20; -15 5 25; 40 -10 20];  % MHz
Exp.Field = 3400;
Exp.mwFreq = 95;
Exp.SampleFrame = [0.3 0.7 1.1];

a = endorfrq(Sys,Exp);
b = endorfrq_perturb(Sys,Exp);

ok = areequal(sort(a),sort(b),1e-4,'rel');
