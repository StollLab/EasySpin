function [err,data] = test(opt,olddata)

% g strain

Sys.S = 3/2;
Sys.D = 300;
Sys.g = [2 2.3];
Sys.lwpp = 2;
Exp.mwFreq = 9.8;
Exp.CenterSweep = [330 160];
Exp.Harmonic = 1;

Opt.Method = 'perturb';
Opt.nKnots = 31;
[x,y] = pepper(Sys,Exp,Opt);

if (opt.Display)
  if isempty(olddata)
    oldy = zeros(size(y));
  else
    oldy = olddata.y;
  end
  Opt.Method = 'matrix';
  [x,y_matrix] = pepper(Sys,Exp,Opt);
  subplot(4,1,[1 2 3]);
  plot(x,y,'r',x,oldy,'g',x,y_matrix,'k');
  legend('pt new','pt old','matrix');
  title(mfilename,'Interpreter','none');
  subplot(4,1,4);
  plot(x,oldy-y,'r',x,y_matrix-y,'k');
  legend('old-new','matrix-new');
end

data.y = y;

if ~isempty(olddata)
  err = ~areequal(y,olddata.y,1);
else
  err = [];
end
