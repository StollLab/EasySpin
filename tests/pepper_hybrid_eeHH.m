function [err,data] = test(opt,olddata)

Sys.S = [3/2 1/2];
Sys.g = [2 2];
Sys.A = [1 0; 0 1]*250;
Sys.Nucs = '1H,1H';
Sys.D = [1 1]*1000;
Sys.ee = 2*1.5*30e3;
Sys.lwpp = 3;

Exp.mwFreq = 9.5;
Exp.Range = [180 500];
Exp.nPoints = 10000;
Exp.CrystalOrientation = [0 pi/4 0];
Exp.Temperature = 40;
Exp.Harmonic = 0;

err = 0;

for S = 1/2:1/2:5/2
  Sys.S(2) = S;
  Opt.Method = 'matrix';
  [x,y1] = pepper(Sys,Exp,Opt);
  Opt.Method = 'hybrid';
  [x,y2] = pepper(Sys,Exp,Opt);
  err = err || ~areequal(y1,y2,0.05*max(y1));
end

data = [];

