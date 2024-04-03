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
[B,spc0] = pepper(Sys,Exp);
Sys.DStrainCorr = +1; % perfect correlation
[B,spcp] = pepper(Sys,Exp);
Sys.DStrainCorr = -1; % anticorrelation
[B,spcm] = pepper(Sys,Exp);

% Plotting
plot(B,spc0,B,spcp,B,spcm);
legend('D/E correlation = 0','D/E correlation = +1','D/E correlation = -1');
legend boxoff
