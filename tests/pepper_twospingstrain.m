function ok = test(opt)

% Comparison of one systems with two very weakly coupled
% g-strain-broadened spins with the sum of their indiviual
% spectra.

% parameters
S1 = 1/2;
S2 = 1/2;
g1 = [2 2.01];
g2 = [2.03 2.04];
gs1 = 4*[1 1 1]*0.0005;
gs2 = 1*[1 1 1]*0.0005;

% combined system
Sys.S = [S1 S2];
Sys.g = [g1; g2];
Sys.ee = 1e-4;
Sys.gStrain = [gs1; gs2];

% system 1
Sys1.S = S1;
Sys1.g = g1;
Sys1.gStrain = gs1;

% system 2
Sys2.S = S2;
Sys2.g = g2;
Sys2.gStrain = gs2;

Exp.mwFreq = 9.5;
Exp.Harmonic = 0;
Exp.CrystalOrientation = rand(1,3)*pi;
Exp.Range = [332 342];

[x,y_both] = pepper(Sys,Exp);
[x,y1] = pepper(Sys1,Exp);
[x,y2] = pepper(Sys2,Exp);
y_both = y_both/2;

if opt.Display
  subplot(2,1,1);
  plot(x,y_both,x,y1+y2);
  legend('both','1+2');
  subplot(2,1,2);
  plot(x,y_both-y1-y2);
  legend('residuals');
end

ok = areequal(y_both,y1+y2,1e-5*max(y_both),'abs');
