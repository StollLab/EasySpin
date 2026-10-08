function ok = test()

% In hybrid mode, nuclei with strains are moved to the exact core, so widths
% equal those from full matrix diagonalization.

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [10 20 30];
Sys.StrainPars = {'A(3)','g(1)'};
Sys.StrainFWHM = [5 0.01];

Exp.mwFreq = 9.5;
Exp.Range = [300 345];
Exp.SampleFrame = [0.2 0.9 0.4];

[B0,~,W0] = resfields(Sys,Exp);
Opt.Hybrid = 1;
Opt.HybridCoreNuclei = [];
[B1,~,W1] = resfields(Sys,Exp,Opt);

ok(1) = areequal(B0,B1,1e-8,'rel');
ok(2) = areequal(W0,W1,1e-8,'rel');
