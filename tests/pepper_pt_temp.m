function [err,data] = test(opt,olddata)

% temperature effects perturbation <-> matrix diagonalization

Sys.S = 1;
Sys.g = [2];
Sys.D = 400;
%Sys.A = [50 80];
%Sys.Nucs = '1H';
Sys.lwpp = 0.5;
Exp.mwFreq = 9.5;
Exp.CenterSweep = [330 60];
Exp.Temperature = 0.5;

Opt.Method = 'matrix';
[x,y0]=pepper(Sys,Exp,Opt);
Opt.Method = 'perturb';
[x,y1]=pepper(Sys,Exp,Opt);

if opt.Display
  plot(x,y0,x,y1);
  legend('matrix','perturb');
  legend boxoff
end
data = [];

err = ~areequal(y0,y1,2e-2*max(y0));
