% Triplet half-field transition and D strain: trimethylenemethane
%===============================================================================
% This script illustrates that the relative intensity of the half-field
% transition and the wings of the allowed transitions in the simulated EPR
% spectrum of a spin triplet is very sensitive to the line broadening model
% included (D strain, HStrain, lwpp).
%
% The spin system is trimethylenemethane (TMM), a non-Kekule diradical with a
% thermally populated triplet ground state, observed in frozen matrices at
% 77 K. Its D3h symmetry gives an axial zero-field splitting, E = 0.
%  Dowd, J. Am. Chem. Soc. 88, 2587 (1966), https://doi.org/10.1021/ja00963a039

clear, clc, clf

D0 = unitconvert(0.024,'cm^-1->MHz');  % center of Gaussian distribution of D values, MHz
Dfwhm = 150;  % full width at half maximum of Gaussian distribution of D values, MHz

Exp.mwFreq = 9.5;  % GHz
Exp.Range = [150 400];  % mT

Opt.GridSize = 90;

Triplet.S = 1;
Triplet.g = 2.0023;
Triplet.D = D0;
Triplet.HStrain = 50;  % MHz, unresolved proton hyperfine couplings

% (1) Simulate spectrum using built-in D strain to model D distribution
TripletStrain = Triplet;
TripletStrain.StrainPars = {'D'};
TripletStrain.StrainFWHM = Dfwhm;
[B,spc_DStrain] = pepper(TripletStrain,Exp,Opt);

% (2) Simulate spectrum using explicit loop over D distribution
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
Triplet.HStrain = 90;  % MHz
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
legend('D strain','loop over D distribution','HStrain only','location','best');
xlabel('magnetic field  (mT)');
ylabel('d\chi''''/dB  (arb.u.)');
