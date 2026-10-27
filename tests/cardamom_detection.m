function ok = test(opt)

% Check harmonic, detection phase and field modulation against chili

Sys.g = [2.009, 2.006, 2.002];
Sys.Nucs = '14N';
Sys.A = unitconvert([6, 36]/10,'mT->MHz'); % MHz
Sys.tcorr = 5e-9; % s
Sys.lw = [0.0, 0.1]; % mT

Par.dtSpatial = 1e-9; % s
Par.dtSpin = Par.dtSpatial;
Par.nSteps = 150;
Par.nTraj = 50;
Par.Model = 'diffusion';

Opt.Method = 'fast';

ExpB.mwFreq = 9.4;
ExpF.Field = 340;
ExpF.mwRange = [9.3 9.8];

Exps = {ExpB,ExpB,ExpB,ExpF};
Exps{1}.Harmonic = 0;
Exps{2}.Harmonic = 0;
Exps{2}.mwPhase = pi/2;
Exps{3}.ModAmp = 0.5;
Exps{4}.mwPhase = pi/2;

for k = numel(Exps):-1:1
  rng(1)
  [x,ycard] = cardamom(Sys,Exps{k},Par,Opt);
  [~,ychili] = chili(Sys,Exps{k});
  ycard = ycard/max(abs(ycard));
  ychili = ychili/max(abs(ychili));
  rmsd = sqrt(mean((ycard-ychili).^2));
  ok(k) = rmsd<0.05;

  if opt.Display
    subplot(numel(Exps),1,k);
    plot(x,ycard,x,ychili);
    legend('cardamom','chili');
    title(sprintf('rmsd %g',rmsd));
  end
end

end
