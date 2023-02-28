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
[x,y1] = pepper(Sys,Exp);

Sys.gFrame = [0 0 0];
Exp.MolFrame = gFrame;
[x,y2] = pepper(Sys,Exp);

if (opt.Display)
  subplot(3,1,[1 2])
  plot(x,y1,x,y2,'r');
  legend('gFrame','MolFrame');
  ylabel('intensity (arb.u.)');
  title('pepper: gFrame vs MolFrame');
  subplot(3,1,3)
  plot(x,y2-y1,'r');
  xlabel('magnetic field (mT)');
end

ok = areequal(y1,y2,1e-3,'rel');
