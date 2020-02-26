function [ok,data] = test(opt,olddata)

% Cu natural abundance mixture

Sys.S = 1/2;
Sys.g = [2 2.2];
Sys.Nucs = 'Cu';
Sys.A = [50 400];
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.Range = [280 350];
[x,y] = pepper(Sys,Exp);

if (opt.Display)
  if ~isempty(olddata)
    subplot(4,1,4);
    plot(x,y-olddata.y);
    subplot(4,1,[1 2 3]);
    plot(x,olddata.y,'k',x,y,'r');
    legend('old','new');
    title('pepper: Cu natural isotope mixture');
    xlabel('magnetic field [mT]');
    ylabel('intensity [a.u.]');
  else
    plot(x,y);
    title('pepper: Cu natural isotope mixture');
    xlabel('magnetic field [mT]');
    ylabel('intensity [a.u.]');
  end
end

data.y = y;

if ~isempty(olddata)
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-5,'abs');
else
  ok = [];
end
