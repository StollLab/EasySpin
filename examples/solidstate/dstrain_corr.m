% Demonstration of correlated D-E strain
%----------------------------------------------------

clear, clf
Sys.S = 1;
Sys.D = [800 80];
Sys.DStrain = [100 33];
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.CenterSweep = [340 100];

% Varying the correlation coefficient
Sys.DStrainCorr = 0; % no correlation
[x,y0] = pepper(Sys,Exp);
Sys.DStrainCorr = +1; % perfect correlation
[x,yp] = pepper(Sys,Exp);
Sys.DStrainCorr = -1; % anticorrelation
[x,ym] = pepper(Sys,Exp);

plot(x,y0,x,yp,x,ym);
legend('correlation = 0','correlation = +1','correlation = -1');
legend boxoff
