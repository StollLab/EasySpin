function ok = test()

% Assert that intensities with a higher-order Zeeman term are the same as
% without, when the term is negligible (electron and nuclear moments combined
% with correct signs)

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [20 60];
SysHo = Sys;
SysHo.Ham312 = [0 0 1e-12 0 0];

Exp.mwFreq = 9.5;
Exp.Range = [300 380];
Exp.MolFrame = [0.3 0.7 0.2];
Exp.CrystalSymmetry = 1;

[~,Int] = resfields(Sys,Exp);
[~,IntHo] = resfields(SysHo,Exp);

ok = areequal(IntHo,Int,1e-6,'rel');
