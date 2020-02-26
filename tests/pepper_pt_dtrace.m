function ok = test(opt)

% compare spectra traceless D and non-traceless D

Sys.S = 1;
Sys.g = 2;
Sys.lwpp = 1;
Exp.mwFreq = 9.8;
Exp.CenterSweep = [350 100];
Exp.Harmonic = 1;

D = [-300 -300 600];
D1 = D + 100;

Sys.D = D;
Opt.Method = 'perturb';
[x1,y1] = pepper(Sys,Exp,Opt);
Sys.D = D1;
Opt.Method = 'perturb';
[x2,y2] = pepper(Sys,Exp,Opt);

if (opt.Display)
  plot(x1,y1,x2,y2);
end

ok = areequal(y1,y2,1e-3,'rel');
