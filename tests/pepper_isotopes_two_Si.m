function [err,data] = test(opt,olddata)

% Si,Si natural abundance mixture, full A tensor specification
Sys.Nucs = 'Si,Si';
Sys.lwpp = 0.1;
Exp.mwFreq = 9.7;
Exp.CenterSweep = [346 10];

A1 = [5 5 5];
A2 = [10 10 10];
Sys.A = [A1; A2];
[x1,y1] = pepper(Sys,Exp);
A1 = diag([5 5 5]);
A2 = diag([10 10 10]);
Sys.A = [A1; A2];
[x2,y2] = pepper(Sys,Exp);

if (opt.Display)
  plot(x1,y1,x2,y2)
  xlabel('magnetic field [mT]');
  ylabel('intensity [a.u.]');
  title('pepper: Si+Si natural abundance mixture');
end

data = [];

err = ~areequal(y1/max(y1),y2/max(y2),1e-5);
