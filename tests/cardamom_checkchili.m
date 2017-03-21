function [err,data] = test(opt,olddata)
% Check results of cardamom against chili

rng(1)

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.tcorr = 50e-9;
Sys.B = 0.34;  % T
Sys.lw = [0.0, 0.1];

Par.dt = 3e-9;
Par.nSteps = ceil(200e-9/Par.dt);
Par.nTraj = 200;

% Exp.mwFreq = mt2mhz(Sys.B, sum(Sys.g)/3);
Exp.mwFreq = 9.4;

% Opt.Verbosity = 1;
Opt.Method = 'Sezer';

[Bcard,ycard] = cardamom(Sys,Par,Exp,Opt);

ycard = ycard/max(ycard);

% Calculate spectrum using chili
% -------------------------------------------------------------------------

[Bchili,ychili] = chili(Sys, Exp);

ychili = ychili/max(ychili);

% Plot for comparison
% -------------------------------------------------------------------------
% figure
% 
% plot(Bcard, ycard, Bchili, ychili)
% ylim([-1.1,1.1])
% ylabel('Im(FFT(M_{+}(t)))')
% xlabel('B (mT)')
% legend('cardamom','chili')
% hold off

rmsd = sqrt(mean((ycard-ychili).^2));

if rmsd < 0.2
  err = 0;
else
  err = 1;
end

data = [];

end
