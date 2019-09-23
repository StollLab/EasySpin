function [err,data] = test(opt,olddata)

% test full gmatrix spec with perturbation theory
Sys.g = [2 0 0; 0 2.05 0; 0 0 2.1];
Sys.lwpp = 1;
Exp.mwFreq = 9.7;
Exp.CenterSweep = [338 40];

Opt.Method = 'matrix';
[x0,y0] = pepper(Sys,Exp,Opt);

Opt.Method = 'perturb1';
[x1,y1] = pepper(Sys,Exp,Opt);

Opt.Method = 'perturb2';
[x2,y2] = pepper(Sys,Exp,Opt);

if opt.Display
  plot(x0,y0,x1,y1,x2,y2);
end

err = ~areequal(y0,y1,5e-3,'abs') || ~areequal(y0,y2,5e-3,'abs');
data = [];
