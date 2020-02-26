function [ok,data] = test(opt,olddata)

% g strain

Sys.g = [2 2.1 2.2];
Sys.lwpp = 2;
Exp.mwFreq = 9.8;
Exp.CenterSweep = [335 60];
Exp.Harmonic = 1;

Opt.Method = 'perturb';
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
  ok = areequal(y/max(y),olddata.y/max(olddata.y),1e-4,'abs');
else
  ok = [];
end
