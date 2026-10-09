% Demonstration of correlated D-E strain
%----------------------------------------------------
% Generic spin triplet (S = 1) with zero-field splitting and strain
% magnitudes typical of organic triplets and biradicals.

clear, clf

Sys.S = 1;
Sys.D = [1500 150];  % D and E, MHz
Sys.StrainPars = {'D(1)','D(2)'};  % D and E
Sys.StrainFWHM = [200 50];  % MHz
Sys.lwpp = 1;

Exp.mwFreq = 9.5;
Exp.CenterSweep = [340 160];

% Varying the correlation coefficient
Sys.StrainCorr = 0; % no correlation
[B,spc0] = pepper(Sys,Exp);
Sys.StrainCorr = +1; % perfect correlation
[B,spcp] = pepper(Sys,Exp);
Sys.StrainCorr = -1; % anticorrelation
[B,spcm] = pepper(Sys,Exp);

% Plotting
plot(B,spc0,B,spcp,B,spcm);
legend('D/E correlation = 0','D/E correlation = +1','D/E correlation = -1');
legend boxoff
