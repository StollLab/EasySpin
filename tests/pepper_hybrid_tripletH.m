function [err,data] = test(opt,olddata)

Sys.S = 1;
Sys.g = 2;
Sys.A = 150;
Sys.Nucs = '1H';
Sys.D = 2000;
Sys.lwpp = 2;

Exp.mwFreq = 9.5;
Exp.Range = [100 440];
Exp.nPoints = 10000;
Exp.Harmonic = 0;

Opt.Method = 'matrix';
[x,yMatrix] = pepper(Sys,Exp,Opt);
Opt.Method = 'hybrid';
[x,yHybrid] = pepper(Sys,Exp,Opt);

if (opt.Display)
  plot(x,yMatrix,'b',x,yHybrid,'r');
  legend('matrix','hybrid');
end

data = [];

err = ~areequal(yMatrix,yHybrid,0.02,'rel');
