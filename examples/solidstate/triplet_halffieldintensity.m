% Triplet half-field transition and D strain
%===============================================================================
% This script illustrates that the relative intensity of the half-field
% transition and the wings of the allowed transitions in the simulated EPR
% spectrum of a spin triplet is very sensitive to the line broadening model
% included (DStrain, HStrain, lwpp).

clear, clc, clf

D0 = 2000;  % center of Gaussian distribution of D values, MHz
Dfwhm = 500;  % full width at half maximum of Gaussian distribution of D values, MHz

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [100 500];  % mT

Opt.GridSize = 90;

Triplet.S = 1;
Triplet.D = D0;
Triplet.DStrain = Dfwhm;
Triplet.HStrain = 50;  % MHz

% (1) Simulate spectrum using build-in DStrain to model D distribution
[B,spc_DStrain] = pepper(Triplet,Exp,Opt);

% (2) Simulate spectrum using explicit loop over D distribution
Triplet.DStrain = 0;
D = linspace(-1,1,51)*2*Dfwhm + D0;
weights = gaussian(D,D0,Dfwhm);
weights = weights/sum(weights);
spc_Dloop = 0;
for k = 1:numel(weights)
  Triplet.D = D(k);
  spc_Dloop = spc_Dloop + weights(k)*pepper(Triplet,Exp,Opt);
end

% (3) Simulate spectrum using HStrain only
Triplet.D = D0;
Triplet.DStrain = 0;
Triplet.HStrain = 280;
[B,spc_HStrain] = pepper(Triplet,Exp,Opt);

% Normalize all spectra
normalize = @(y)y/max(y);
spc_DStrain = normalize(spc_DStrain);
spc_Dloop = normalize(spc_Dloop);
spc_HStrain = normalize(spc_HStrain);

% Plotting
plot(B,spc_DStrain,B,spc_Dloop,B,spc_HStrain);
grid on
axis tight
legend('DStrain','loop over D distribution','HStrain only','location','best');
xlabel('magnetic field  (mT)');
ylabel('d\chi''''/dB  (arb.u.)');
