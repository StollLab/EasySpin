function [ok,data] = test(opt,olddata)

% g strain

Sys = struct('S',1/2,'g',[2 2.1 2.2],'gStrain',[1 2 3]*0.01);
Exp = struct('mwFreq',9.5,'Range',[290 350]);
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
  title('pepper: g strain');
end

data.y = y;

if ~isempty(olddata)
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-4,'abs');
else
  ok = [];
end
