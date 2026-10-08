function ok = test()

% Opt.Freq2Field=0 omits the 1/g factor from intensities and returns widths
% in MHz instead of mT

Sys.g = [2 2.1 2.2];
Sys.Nucs = '1H';
Sys.A = [10 20 30];
Sys.HStrain = [5 10 15];
Sys.gStrain = [0.01 0.02 0.03];
Sys.AStrain = [2 3 4];
Exp.mwFreq = 9.5;
Exp.Range = [280 360];
Exp.SampleFrame = [0 0 0; 0 pi/2 0];  % field along z(Mol) and x(Mol)

[P1,I1,W1] = resfields_perturb(Sys,Exp);
Opt.Freq2Field = 0;
[P0,I0,W0] = resfields_perturb(Sys,Exp,Opt);

% 1/g factor for each orientation, in mT/MHz
geff = Sys.g([3 1]);
dBdE = planck/bmagn*1e9./geff;

ok(1) = areequal(P0,P1,1e-10,'abs');
ok(2) = areequal(I0.*dBdE,I1,1e-10,'rel');
ok(3) = areequal(W0.*dBdE,W1,1e-10,'rel');
