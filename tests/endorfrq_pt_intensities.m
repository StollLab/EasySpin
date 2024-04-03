function ok = test()

% Compare intensites between matrix diagonalization and
% perturbation theory

OptMD.Method='matrix';
OptMD.separate='transitions';
OptMD.Threshold=0.3;

OptST.Method='perturb2';
OptST.separate='transitions';
OptST.PerturbOrder=2;

MnBsp.g = 1.96;
MnBsp.Nucs = '63Cu';
MnBsp.A = [250];
MnBsp.HStrain=10;

Endor.Field=1234.5;
Endor.Range=[100 150];
Endor.SampleFrame = [0 0 0];

[PosMD,IntMD]=endorfrq(MnBsp,Endor,OptMD);
[PosMD,idx] = sort(PosMD); IntMD = IntMD(idx);
[PosST,IntST]=endorfrq_perturb(MnBsp,Endor,OptST);
[PosST,idx] = sort(PosST); IntST = IntST(idx);

ok = areequal(IntMD,IntST,1e-3,'rel');
