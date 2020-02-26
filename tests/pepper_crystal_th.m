function [ok,data] = test(opt,olddata)

% Crystal with Th symmetry and rhombic g

clear Sys Exp Opt
Sys.g = [2.0 2.1 2.2];
Sys.gFrame = [-78 -40 -30]*pi/180;
Sys.lw = 0.5;
Exp.mwFreq = 9.5;
Exp.Range = [290 350];

n = [1;2;3];
[phi, theta] = vec2ang(n);
Exp.CrystalOrientation = [phi theta 0];
Exp.CrystalSymmetry = 'Th';
Opt.Verbosity = 0;

[x,y] = pepper(Sys,Exp,Opt);

if (opt.Display)
  if ~isempty(olddata)
    subplot(4,1,[1 2 3]);
    plot(x,olddata.y,'k',x,y,'r');
    legend('old','new');
    legend boxoff
    subplot(4,1,4);
    plot(x,y-olddata.y);
  else
    plot(x,y);
  end
  xlabel('magnetic field (mT)');
  ylabel('intensity (arb.u.)');
  title('pepper: Th crystal');
end

data.y = y;

if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-4,'abs');
else
  ok = [];
end
