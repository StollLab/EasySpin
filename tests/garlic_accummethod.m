function ok = test(opt)

% Test whether all accumulation methods in Opt.AccumMethod yield the
% same result.

Sys.Nucs = '1H';
Sys.A = 10;
Sys.lwpp = [0 0.01];
Exp.mwFreq = 9.5;
Exp.Range = [337.5 340.5];
Exp.Harmonic = 0;
Exp.nPoints = 1e4;

Opt.AccumMethod = 'explicit';
[x,y0] = garlic(Sys,Exp,Opt);
Opt.AccumMethod = 'template';
[x,y1] = garlic(Sys,Exp,Opt);
Opt.AccumMethod = 'binning';
[x,y2] = garlic(Sys,Exp,Opt);

if opt.Display
  plot(x,y0,x,y1,x,y2);
  legend('explicit','template','binning');
  legend boxoff
end

thr = 2e-2;
ok = areequal(y0,y1,thr,'rel') && areequal(y0,y2,thr,'rel');
