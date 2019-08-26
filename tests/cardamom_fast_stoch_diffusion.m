function [err,data] = test(opt,olddata)
% Check results of cardamom against chili

rng(1)

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10); % MHz
Sys.tcorr = 5e-9; % s
Sys.lw = [0.0, 0.1]; % mT

Par.dt = 1e-9; % s
Par.Dt = Par.dt;
Par.nSteps = 150;
Par.nTraj = 50;
Par.Model = 'diffusion';

Exp.mwFreq = 9.4;

Opt.Method = 'fast';

% Calculate spectrum using cardamom
[Bcard,ycard] = cardamom(Sys,Exp,Par,Opt);
ycard = ycard/max(ycard);

% Calculate spectrum using chili
[Bchili,ychili] = chili(Sys, Exp);

ychili = ychili/max(ychili);

rmsd = sqrt(mean((ycard-ychili).^2));
err = rmsd>0.05;

if opt.Display
  plot(Bcard, ycard, Bchili, ychili)
  ylim([-1.1,1.1])
  ylabel('Im(FFT(M_{+}(t)))')
  xlabel('B (mT)')
  legend('cardamom','chili')
  hold off
end

data = [];

end
