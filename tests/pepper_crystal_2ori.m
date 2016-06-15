function [err,data] = test(opt,olddata)

% Crystal with D2h symmetry and rhombic g
% several orientations

Sys.g = [2 2.1 2.2];
Sys.gFrame = -[78 40 30]*pi/180;
Sys.lw = 0.5;
Exp.mwFreq = 9.5;
Exp.Range = [290 350];

Exp.CrystalOrientation = vec2ang([1 2 3; 1 -2 4; 0 0 1; 5 2 -3]).';
Exp.CrystalOrientation(:,3) = 0;
Exp.CrystalSymmetry = 'D2h';

[x,y] = pepper(Sys,Exp);

if (opt.Display)
  if ~isempty(olddata)
    subplot(4,1,[1 2 3]);
    h = plot(x,olddata.y,'k',x,y,'r');
    legend('old','new');
    %set(h(1),'LineWidth',2);
    subplot(4,1,4);
    plot(x,y-olddata.y);
  else
    plot(x,y);
  end
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Th crystal');
end

data.y = y;

if ~isempty(olddata)
  err = ~areequal(y/max(y),olddata.y/max(olddata.y),1e-4);
else
  err = [];
end
