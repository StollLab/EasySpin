function ok = test(opt)

% Check that the template accumulation gives the same dispersion spectrum as
% the explicit evaluation of the line shapes

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 20;
Sys.lw = [0 0.1];
Exp.mwFreq = 9.5;
Exp.Range = [336 343];
Exp.nPoints = 2000;
Exp.mwPhase = pi/2;

harmonics = [0 1];
for k = numel(harmonics):-1:1
  Exp.Harmonic = harmonics(k);
  Opt.AccumMethod = 'explicit';
  [x,y0] = garlic(Sys,Exp,Opt);
  Opt.AccumMethod = 'template';
  [x,y1] = garlic(Sys,Exp,Opt);
  ok(k) = areequal(y0,y1,1e-2,'rel');

  if opt.Display
    subplot(numel(harmonics),1,k);
    plot(x,y0,x,y1,'--');
    legend('explicit','template');
  end
end

end
