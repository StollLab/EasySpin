function [err,data] = test(opt,olddata)
% Check results of cardamom against chili

rng(1)

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = [6, 36]/1e4;
Sys.tcorr = 50e-9;
Sys.B = 0.34;  % T
Sys.lw = [0.0, 0.1];

Par.dt = 3e-9;
Par.nSteps = ceil(200e-9/Par.dt);
Par.nTraj = 200;

Exp.mwFreq = mt2mhz(Sys.B, sum(Sys.g)/3);
% Exp.mwFreq = 9.4;  % GHz

Opt.Verbosity = 1;
Opt.Method = 'DeSensi';

[freq,spc,expval] = cardamom(Sys,Par,Exp,Opt);

xcard = freq/1e6;
ycard = spc/max(spc);

% Calculate spectrum using chili
% -------------------------------------------------------------------------
Sys.A = mt2mhz([6, 36]/10);
Sys.lw = [0.0, 0.1];
% 
[B,spc] = chili(Sys, Exp);
% 
xchili = mt2mhz(B-(max(B)+min(B))/2);
ychili = spc/max(spc);
% 
% filter = xcard<max(xchili) & xcard>min(xchili);
% xcard = xcard(filter);
% ycard = ycard(filter);
% 
ycard = interp1(xcard,ycard,xchili);

% Plot for comparison
% -------------------------------------------------------------------------
figure

plot(B, ycard, B, ychili)
ylim([-1.1,1.1])
xlim([min(B),max(B)])
ylabel('Im(FFT(M_{+}(t)))')
xlabel('f (MHz)')
legend('cardamom','chili')
hold off

rmsd = sqrt(mean((ycard-ychili).^2));

if rmsd < 0.2
  err = 0;
else
  err = 1;
end

data = [];

end
