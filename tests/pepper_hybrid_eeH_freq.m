function [err,data] = test(opt,olddata)

% Frequency sweeps: Hybrid vs. full matrix diagonalization
%--------------------------------------------------------------

Sys.S = [3/2 2];
Sys.g = [2 2];
Sys.A = [1 0]*100;
Sys.Nucs = '1H';
Sys.D = [1 1]*1000;
Sys.ee = 2*1.5*30e3;
Sys.lwpp = 20;

Exp.Field = 350;
Exp.mwRange = [7 11];
Exp.nPoints = 40000;
Exp.CrystalOrientation = [0 0 0] ;
Exp.Temperature = 8;

Opt.Method = 'matrix';
[x1,y1] = pepper(Sys,Exp,Opt);
Opt.Method = 'hybrid';
[x2,y2] = pepper(Sys,Exp,Opt);

if (opt.Display)
  plot(x1,y1,'b',x2,y2,'r');
  legend('matrix','hybrid');
end

data = [];

err = ~areequal(y1,y2,0.02*max(y1));
