function [err,data] = test(opt,olddata)

% For a high-spin system, compare line positions and intensities
% as obtained by matrix diagonalization and by perturbation theory.

clear Sys Exp
Sys.S = 5/2;
Sys.D = 10;
Sys.lwpp = 0.3;

Exp.mwFreq = 9.5;
Exp.Range = [300 380];
Exp.CrystalOrientation = [pi/4 pi/4 0];

[Ba,Ia] = resfields_perturb(Sys,Exp);
[Bb,Ib] = resfields(Sys,Exp);

[Ba,idx] = sort(Ba); Ia = Ia(idx);
[Bb,idx] = sort(Bb); Ib = Ib(idx);

err(1) = any(abs(Ba-Bb)>0.0001);
err(2) = any(abs(Ia-Ib)>0.05);

data = [];
