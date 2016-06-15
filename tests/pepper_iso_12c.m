function [err,data] = test(opt,olddata)

% 12C natural abundance mixture

Sys.S = 1/2;
Sys.g = [2 2.2];
Sys.Nucs = '12C';
Sys.A = [50 400];
Sys.lwpp = 0.3;

Exp.mwFreq = 9.5;
Exp.Range = [280 350];

[x,y] = pepper(Sys,Exp);

if (opt.Display)
  if ~isempty(olddata)
    subplot(4,1,[1 2 3]);
    h = plot(x,olddata.y,'k',x,y,'r');
    legend('old','new');
    subplot(4,1,4);
    plot(x,y-olddata.y);
  else
    plot(x,y);
  end
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Cu+N natural isotope mixture');
end

data.y = y;

if ~isempty(olddata)
  err = ~areequal(y/max(y),olddata.y/max(olddata.y),1e-5);
else
  err = [];
end
