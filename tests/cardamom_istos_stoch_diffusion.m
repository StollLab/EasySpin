function ok = test(opt)

% Check results of cardamom against chili

rng(1)

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.tcorr = 100e-9;
Sys.B = 0.34;  % T
Sys.lw = [0.0, 0.1];

Par.dtSpatial = 3e-9;
Par.dtSpin = 3e-9;
Par.nSteps = ceil(100e-9/Par.dtSpin);
Par.nTraj = 50;
Par.Model = 'diffusion';

Exp.mwFreq = 9.4;

Opt.Verbosity = 0;
Opt.Method = 'ISTOs';

[Bcard,ycard] = cardamom(Sys,Exp,Par,Opt);

ycard = ycard/max(ycard);

% Calculate spectrum using chili
% -------------------------------------------------------------------------

[Bchili,ychili] = chili(Sys, Exp);

ychili = ychili/max(ychili);

% Plot for comparison
% -------------------------------------------------------------------------

rmsd = sqrt(mean((ycard-ychili).^2));

ok = rmsd<0.2;

if opt.Display
  plot(Bcard, ycard, Bchili, ychili)
  ylim([-1.1,1.1])
  ylabel('Im(FFT(M_{+}(t)))')
  xlabel('B (mT)')
  legend('cardamom','chili')
  hold off
end

end
