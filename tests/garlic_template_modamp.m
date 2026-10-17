function ok = test(opt)

% Test whether template accumulation gives the same result as explicit
% accumulation when field modulation (Exp.ModAmp) is used.

Sys.Nucs = '1H';
Sys.A = 10;
Sys.lwpp = [0 0.05];
Exp.mwFreq = 9.5;
Exp.Range = [337.5 340.5];
Exp.ModAmp = 0.1;
Exp.nPoints = 4000;

for h = 1:2
  Exp.Harmonic = h;
  Opt.AccumMethod = 'explicit';
  [x,y0] = garlic(Sys,Exp,Opt);
  Opt.AccumMethod = 'template';
  [x,y1] = garlic(Sys,Exp,Opt);

  if opt.Display
    subplot(2,1,h);
    plot(x,y0,x,y1);
    legend('explicit','template');
    legend boxoff
  end

  ok(h) = areequal(y0,y1,2e-2,'rel');
end
