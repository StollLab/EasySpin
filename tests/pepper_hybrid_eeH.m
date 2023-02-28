function ok = test(opt)

% g strain

Sys.S = [3/2 2];
Sys.g = [2 2];
Sys.A = [1 0]*100;
Sys.Nucs = '1H';
Sys.D = [1 1]*1000;
Sys.ee = 2*1.5*30e3;
Sys.lwpp = 2;

Exp.mwFreq = 9.5;
Exp.Range = [180 500];
Exp.nPoints = 10000;
Exp.SampleFrame = [0 0 0];
Exp.Temperature = 8;

Opt.Method = 'matrix';
[x,y1] = pepper(Sys,Exp,Opt);
Opt.Method = 'hybrid';
[x,y2] = pepper(Sys,Exp,Opt);

if (opt.Display)
  plot(x,y1,'b',x,y2,'r');
  legend('matrix','hybrid');
end

ok = areequal(y1,y2,0.02,'rel');
