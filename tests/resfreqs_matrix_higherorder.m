function ok = test()

% Assert that resfreqs_matrix works with a higher-order Zeeman term and gives
% the same result as without, when the term is negligible

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [20 60];
SysHo = Sys;
SysHo.Ham312 = [0 0 1e-12 0 0];

Exp.Field = 340;
Exp.mwRange = [9 10];
Exp.MolFrame = [0.3 0.7 0.2];
Exp.CrystalSymmetry = 1;

[nu,Int] = resfreqs_matrix(Sys,Exp);
[nuHo,IntHo] = resfreqs_matrix(SysHo,Exp);

ok = areequal(nuHo,nu,1e-6,'rel') && areequal(IntHo,Int,1e-6,'rel');
