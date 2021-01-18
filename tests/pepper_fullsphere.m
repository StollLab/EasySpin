function ok = test(opt)

% Simulations using D2h vs Ci vs C1 grid symmetries

Sys.g = [2 2.05 2.2];
Sys.lw = 1;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.Range = [305 345];

Opt.GridSize = [10 2];

Opt.GridSymmetry = 'D2h';
[B,spc0] = pepper(Sys,Exp,Opt);
Opt.GridSymmetry = 'Ci';
[B,spc1] = pepper(Sys,Exp,Opt);
Opt.GridSymmetry = 'C1';
[B,spc2] = pepper(Sys,Exp,Opt);

scale = max(spc0);
spc0 = spc0/scale;
spc1 = spc1/scale;
spc2 = spc2/scale;

threshold = 0.03;
ok(1) = areequal(spc1,spc0,threshold,'abs');
ok(2) = areequal(spc2,spc0,threshold,'abs');

if opt.Display
  plot(B,spc0,B,spc1,B,spc2);
  legend('D2h','Ci','C1');
end
