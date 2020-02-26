function ok = test(opt)

Sys.S = 3/2;
Sys.D = 300;


Sys.g = [2 2.2];
Sys.lw = 1;
%Sys.Nucs = '1H';
%Sys.A = [500 400];
Exp.mwFreq = 9.5;
Exp.Range = [250 380];
Exp.Harmonic = 0;
Opt.nKnots=91;

Opt.Method = 'matrix';
[x0,y0]=pepper(Sys,Exp,Opt);
Opt.Method = 'perturb2';
[x1,y1]=pepper(Sys,Exp,Opt);


if opt.Display
  plot(x0,y0,x1,y1);
  legend('matrix','perturb2');
end

ok = areequal(y0,y1,0.01,'rel');
