function ok = test(opt)

% Compare approximate and explicit field-swept spectra (general)

Sys.g = [2.01 2.003];
Sys.lw = 0.1;
Sys.tcorr = 1e-9;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.nPoints = 200;
Exp.Range = [337 339];

Opt.LiouvMethod = 'general';

Opt.FieldSweepMethod = 'approxlin';
[x,y1] = chili(Sys,Exp,Opt);

Opt.FieldSweepMethod = 'explicit';
[x,y2] = chili(Sys,Exp,Opt);

if opt.Display
  plot(x,y1,x,y2);
  legend('approximate','explicit');
end

ok = areequal(y1/max(y1),y2/max(y2),0.01,'abs');
