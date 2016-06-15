function [err,data] = test(opt,olddata)

% g strain

Sys.S = 1;
Sys.g = 2;
Sys.D = 600;
Sys.lwpp = 1;
Exp.mwFreq = 9.8;
Exp.CenterSweep = [350 100];
Exp.Harmonic = 1;

Opt.Method = 'perturb';
[x,y] = pepper(Sys,Exp,Opt);
y = y/max(y);

if (opt.Display)
  if isempty(olddata)
    oldy = zeros(size(y));
  else
    oldy = olddata.y/max(olddata.y);
  end
  Opt.Method = 'matrix';
  [x,y2] = pepper(Sys,Exp,Opt);
  y2 = y2/max(y2);
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
  err = ~areequal(y/max(y),olddata.y/max(olddata.y),1e-4);
else
  err = [];
end
