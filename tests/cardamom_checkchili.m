function [err,data] = test(opt,olddata)
% Check results of cardamom against chili

%% Calculate spectrum using cardamom
Sys.Nucs = '14N';

Sys.g = [2.00809, 2.00585, 2.00202];
Sys.A = [6.2, 4.3, 36.9]/1e4;
Sys.tcorr = 50e-9;
% Sys.B = 0.34;  % T
Sys.lw = [0.0, 0.01];

Par.dt = 3e-9;
Par.nSteps = ceil(200e-9/Par.dt);
Par.nTraj = 200;

% Exp.mwFreq = mt2mhz(Sys.B, sum(Sys.g)/3);
Exp.mwFreq = 9.4;  % GHz

Opt.Verbosity = 1;

[freq,spc] = cardamom(Sys,Par,Exp,Opt);

xcard = freq/1e6;
ycard = spc/max(spc);

%% Calculate spectrum using chili
Sys.A = mt2mhz([6.2, 4.3, 36.9]/10);
% Sys.lw = [0.0, 0.01];

[B,spc] = chili(Sys, Exp);

xchili = mt2mhz(B-(max(B)+min(B))/2);
ychili = spc/max(spc);

% filter = xnum<max(xchili) & xnum>min(xchili);
% xnum = xnum(filter);
% ynum = ynum(filter);

ycard = interp1(xcard,ycard,xchili);

%% Plot for comparison
figure

% subplot(2,1,1);
% % plot(t/1e-6, imag(expval_avg))
% plot(t/1e-6, imag(sum(cell2mat(expval),2)))
% ylabel('Im(M_{+}(t))')
% xlabel('t ({\mu}s)')

% subplot(2,1,2);
% plot(mhz2mt(xcard)+340, ycard, B, ychili)
plot(B, ycard, B, ychili)
ylim([-1.1,1.1])
xlim([min(B),max(B)])
ylabel('Im(FFT(M_{+}(t)))')
xlabel('f (MHz)')
legend('cardamom','chili')

rmsd = sqrt(mean((ycard-ychili).^2));

if rmsd < 0.2
  err = 0;
else
  err = 1;
end

data = [];

end
