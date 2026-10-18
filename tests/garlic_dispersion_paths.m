function ok = test(opt)

% Check that explicit and binning accumulation give the same dispersion
% spectrum when both Gaussian and Lorentzian broadening are present. This
% catches a phase that is applied twice.

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 20;
Sys.lw = [0.1 0.05];
Exp.mwFreq = 9.5;
Exp.Range = [336 343];
Exp.nPoints = 2000;

harmonics = [0 1];
for k = numel(harmonics):-1:1
  Exp.Harmonic = harmonics(k);
  Exp.mwPhase = pi/2;
  Opt.AccumMethod = 'explicit';
  [x,y0] = garlic(Sys,Exp,Opt);
  Opt.AccumMethod = 'binning';
  [x,y1] = garlic(Sys,Exp,Opt);
  % Binning rounds the line positions to the grid, which gives deviations
  % of a few percent. A phase applied twice would invert the spectrum.
  idx = round(0.05*numel(x)):round(0.95*numel(x));
  ok(k) = areequal(y0(idx),y1(idx),5e-2,'rel');

  if opt.Display
    subplot(numel(harmonics),1,k);
    plot(x,y0,x,y1,'--');
    legend('explicit','binning');
  end
end

end
