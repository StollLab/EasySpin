function [err,data] = test(opt,olddata)

% Cu natural abundance mixture

Sys.S = 1/2;
Sys.g = [2 2.2];
Sys.Nucs = 'Cu,N';
Sys.A = [50 400; 30 30];
Sys.lwpp = 0.3;

Exp.mwFreq = 9.5;
Exp.Range = [280 350];
Exp.nPoints = 1e4;

Opt.nKnots = 91;
Opt.Method = 'perturb';
[x,y] = pepper(Sys,Exp,Opt);

if (opt.Display)
  if ~isempty(olddata)
    subplot(4,1,[1 2 3]);
    y0 = olddata.y/max(olddata.y);
    y1 = rescale(y,y0,'lsq0');
    h = plot(x,y0,'k',x,y1,'r');
    legend('old','new');
    subplot(4,1,4);
    plot(x,y1-y0);
  else
    plot(x,y);
  end
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Cu+N natural isotope mixture');
end

data.y = y;

if ~isempty(olddata)
  err = ~areequal(y,olddata.y,1e-4);
else
  err = [];
end
