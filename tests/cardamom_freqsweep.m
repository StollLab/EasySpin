function ok = test(opt)

% Check results of cardamom frequency sweep against chili

rng(1)

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.tcorr = 5e-9;
Sys.B = 0.34;  % T
Sys.lw = [0.0, 0.1];

Par.dtSpatial = 1e-9;
Par.dtSpin = 1e-9;
Par.nSteps = ceil(150e-9/Par.dtSpin);
Par.nTraj = 50;
Par.Model = 'diffusion';

Exp.Field = 340;
Exp.mwRange = [9.4, 9.8];

Opt.Verbosity = 0;
Opt.Method = 'fast';
Opt.speccon = 1;

[fcard,ycard] = cardamom(Sys,Exp,Par,Opt);

ycard = ycard/max(ycard);

% Calculate spectrum using chili
% -------------------------------------------------------------------------

[fchili,ychili] = chili(Sys, Exp);

ychili = ychili/max(ychili);

% Plot for comparison
% -------------------------------------------------------------------------

rmsd = sqrt(mean((ycard-ychili).^2))/(max(ychili)-min(ychili));

ok = rmsd<0.1;

if opt.Display
  plot(fcard, ycard, fchili, ychili)
  ylim([-1.1,1.1])
  ylabel('Im(FFT(M_{+}(t)))')
  xlabel('f (GHz)')
  legend('cardamom','chili')
  hold off
end

end
