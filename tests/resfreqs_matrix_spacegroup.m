function ok = test(opt)

% Make sure frequency-swept crystal EPR simulation does not crash
% in the presence of a space group with more than one site.

Sys.g = [2.008 2.006 2.002];
Sys.lwpp = 1;
Sys.gFrame = [40 30 20]*pi/180;

Exp1.Field = 1100;
Exp1.mwRange = [30.7 31];
Exp1.CrystalSymmetry = 220;
Exp1.MolFrame = [10 20 30]*pi/180;
Exp1.SampleFrame = [40 70 20]*pi/180;

[nu,spc] = pepper(Sys,Exp1);

if opt.Display
  plot(nu,spc);
end

ok = true;
