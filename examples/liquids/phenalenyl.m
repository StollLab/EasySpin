% Phenalenyl radical anion, including 13C satellite lines
%==========================================================================
% EasySpin can automatically simulate the spectra resulting from
% isotope mixtures. Here, we simulate the spectra of the phenalenyl
% radicals with (a) pure 12C and (b) natural-abundance 13C. The small
% 13C content leads to small satellite peaks.

% Hyperfine values taken from Gerson book, Table 8.4, page 225

clear, clc

% Spin system parameters
%--------------------------------------------------------------------------
AH1 = -0.629;  % in mT, (H1,H3,H4,H6,H7,H9)
AH2 = +0.181;  % in mT, (H2,H5,H8)
AC1 = +0.966;  % in mT, (C1,C3,C4,C6,C7,C9)
AC2 = -0.784;  % in mT, (C2,C5,C8,C3a,C6a,C9a)
Sys.g = 2;
Sys.n = [6 3 6 6];
Sys.A = unitconvert([AH1 AH2 AC1 AC2],'mT->MHz');  % mT -> MHz
Sys.lwpp = [0, 0.01];  % Lorentzian line shape

% Experimental parameters
%--------------------------------------------------------------------------
Exp.mwFreq = 9.5;  % GHz
Exp.CenterSweep = [339.4 8];  % mT
Exp.nPoints = 10000;

% Simulations
%--------------------------------------------------------------------------
% Spectrum of species with all 12C
Sys.Nucs = '1H,1H,12C,12C';
[B,spc0] = garlic(Sys,Exp);

% Spectra with species containing natural-abundance 13C
Sys.Nucs = '1H,1H,C,C';
[B,spc1] = garlic(Sys,Exp);

% Plotting
%--------------------------------------------------------------------------
plot(B,spc0,B,spc1);
xlabel('magnetic field (mT)');
axis tight
legend('12C only','12C+13C');
legend boxoff
