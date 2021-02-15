function [ok,data] = test(opt,olddata)

% g strain

Sys.S = 1/2;
Sys.g = [2 2.2];
Sys.Nucs = '63Cu,1H';
Sys.A = [70 500; 50 50];
Sys.lwpp = 0.4;
Exp.mwFreq = 9.8;
Exp.CenterSweep = [330 100];
Exp.Harmonic = 1;

Opt.Method = 'perturb';
Opt.GridSize = 31;
[x,y] = pepper(Sys,Exp,Opt);

if (opt.Display)
  if isempty(olddata)
    oldy = zeros(size(y));
  else
    oldy = olddata.y;
  end
  Opt.Method = 'matrix';
  [x,y2] = pepper(Sys,Exp,Opt);
  subplot(4,1,[1 2 3]);
  plot(x,y,'r',x,oldy,'g',x,y2,'k');
  legend('pt new','pt old','matrix');
  title(mfilename,'Interpreter','none');
  subplot(4,1,4);
  plot(x,oldy-y,'r',x,y2-y,'k');
  legend('old-new','matrix-new');
end

data.y = y;

if ~isempty(olddata)
  ok = areequal(y,olddata.y,1e-4,'abs');
else
  ok = [];
end
