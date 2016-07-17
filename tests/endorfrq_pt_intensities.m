function [err,data] = test(opt,olddata)

%-------------------------------------------------------------
% Compare intensites between matrix diagonalization and
% perturbation theory
%-------------------------------------------------------------

OptMD.Method='matrix';
OptMD.Output='separate';
OptMD.Threshold=0.3;

OptST.Method='perturb2';
OptST.Output='separate';
OptST.PerturbOrder=2;

MnBsp.g = 1.96;
MnBsp.Nucs = '63Cu';
MnBsp.A = [250];
MnBsp.HStrain=10;

Endor.Field=1234.5;
Endor.Range=[100 150];
Endor.CrystalOrientation = [0 0 0];

[PosMD,IntMD]=endorfrq(MnBsp,Endor,OptMD);
[PosMD,idx] = sort(PosMD); IntMD = IntMD(idx);
[PosST,IntST]=endorfrq_perturb(MnBsp,Endor,OptST);
[PosST,idx] = sort(PosST); IntST = IntST(idx);

err = any(abs(IntMD./IntST)-1>1e-3);

data = [];
