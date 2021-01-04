function ok = test(opt)

% Check results of cardamom against chili

rng(1)

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.g = [2.009, 2.006, 2.002];
Sys.Nucs = '14N';
Sys.A = mt2mhz([6, 36]/10); % MHz

Sys.tcorr = 5e-9; % s
Sys.lw = [0.0, 0.1]; % mT

Exp.mwFreq = 9.4;

Par.dtSpatial = 1e-9; % s
Par.dtSpin = Par.dtSpatial;
Par.nSteps = 150;
Par.nTraj = 50;
Par.Model = 'diffusion';

Opt.Method = 'fast';

% Calculate spectrum using cardamom
[Bcard,ycard] = cardamom(Sys,Exp,Par,Opt);
ycard = ycard/max(ycard);

% Calculate spectrum using chili
[Bchili,ychili] = chili(Sys, Exp);

ychili = ychili/max(ychili);

rmsd = sqrt(mean((ycard-ychili).^2));
ok = rmsd<=0.05;

if opt.Display
  plot(Bcard, ycard, Bchili, ychili)
  ylim([-1.1,1.1])
  ylabel('Im(FFT(M_{+}(t)))')
  xlabel('B (mT)')
  legend('cardamom','chili')
  hold off
end

end
