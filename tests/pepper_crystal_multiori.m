function [ok,data] = test(opt,olddata)

% Crystal with D2h symmetry and rhombic g, several orientations

Sys.g = [2 2.1 2.2];
Sys.gFrame = -[78 40 30]*pi/180;
Sys.lw = 0.5;
Exp.mwFreq = 9.5;
Exp.Range = [290 350];

Exp.CrystalSymmetry = 'D2h';
[phi,theta] = vec2ang([1 2 3; 1 -2 4; 0 0 1; 5 2 -3]);
chi = zeros(size(phi(:)));
Exp.SampleFrame = [-chi -theta(:) -phi(:)];

[B,spc] = pepper(Sys,Exp);

if opt.Display
  if ~isempty(olddata)
    subplot(4,1,[1 2 3]);
    plot(B,olddata.y,B,spc);
    legend('old','new');
    subplot(4,1,4);
    plot(B,spc-olddata.y);
  else
    plot(B,spc);
  end
  xlabel('magnetic field (mT)');
  ylabel('intensity (arb.u.)');
end

data.y = spc;

if ~isempty(olddata)
  ok = areequal(spc,olddata.y,1e-4,'rel');
else
  ok = [];
end
