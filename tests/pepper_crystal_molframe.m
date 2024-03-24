function ok = test(opt)

% Crystal with Th symmetry and rhombic g

clear Sys Exp Opt
g = [2.0 2.1 2.2];
gFrame = [30 40 78]*pi/180;

Sys.g = g;
Sys.lw = 0.5;
Exp.mwFreq = 9.5;
Exp.Range = [290 350];

Exp.SampleFrame = rand(1,3)*pi;
Exp.CrystalSymmetry = 130;

Sys.gFrame = gFrame;
Exp.MolFrame = [0 0 0];
[B,spc_gFrame] = pepper(Sys,Exp);

Sys.gFrame = [0 0 0];
Exp.MolFrame = gFrame;
[B,spc_MolFrame] = pepper(Sys,Exp);

if opt.Display
  subplot(3,1,[1 2])
  plot(B,spc_gFrame,B,spc_MolFrame,'r');
  legend('gFrame','MolFrame');
  ylabel('intensity (arb.u.)');
  title('pepper: gFrame vs. MolFrame');
  subplot(3,1,3)
  plot(B,spc_MolFrame-spc_gFrame,'r');
  xlabel('magnetic field (mT)');
end

ok = areequal(spc_gFrame,spc_MolFrame,1e-3,'rel');
