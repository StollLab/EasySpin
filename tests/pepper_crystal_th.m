function [ok,data] = test(opt,olddata)

% Crystal with Th symmetry and rhombic g

Sys.g = [2.0 2.1 2.2];
Sys.gFrame = [-78 -40 -30]*pi/180;
Sys.lw = 0.5;
Exp.mwFreq = 9.5;
Exp.Range = [290 350];

n = [1; 2; 3];
[phi, theta] = vec2ang(n);
Exp.SampleFrame = [0 -theta -phi];
Exp.CrystalSymmetry = 'Th';
Opt.Verbosity = 0;

[B,spc] = pepper(Sys,Exp,Opt);

if opt.Display
  if ~isempty(olddata)
    subplot(4,1,[1 2 3]);
    plot(B,olddata.y,'k',B,spc,'r');
    legend('old','new');
    legend boxoff
    subplot(4,1,4);
    plot(B,spc-olddata.y);
  else
    plot(B,spc);
  end
  xlabel('magnetic field (mT)');
  ylabel('intensity (arb.u.)');
  title('pepper: Th crystal');
end

data.y = spc;

if ~isempty(olddata)
  ok = areequal(spc,olddata.y,1e-5,'rel');
else
  ok = [];
end
