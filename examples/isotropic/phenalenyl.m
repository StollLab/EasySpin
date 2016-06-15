% Phenalenyl radical anion, including 13C satellite lines
%==========================================================================
% EasySpin can automatically simulate the spectra resulting from
% isotope mixtures. Here, we simulate the spectra of the phenalenyl
% radicals with (a) pure 12C and (b) natural-abundance 13C. The small
% 13C contents leads to small satellite peaks.

% Hyperfine values taken from Gerson book, Table 8.4, page 225

clear, clf

% Parameters
%--------------------------------------------------------------------------
A_H = [-0.629 +0.181]; % in mT
A_C = [+0.966 -0.784 -0.784]; % in mT
Sys.g = 2;
Sys.n = [6 3 1 1 1];
Sys.A = mt2mhz([A_H A_C]);
Sys.lwpp = [0, 0.01];  % Lorentzian line shape

Exp.mwFreq = 9.5;
Exp.CenterSweep = [339.4 8];
Exp.nPoints = 10000;

% Simulations
%--------------------------------------------------------------------------
% Spectrum of species with all 12C 
Sys.Nucs = '1H,1H,12C,12C,12C';
[x,y0] = garlic(Sys,Exp);

% Spectra with species containing natural-abundance 13C
Sys.Nucs = '1H,1H,C,C,C';
[x,y1] = garlic(Sys,Exp);

% Plotting
%--------------------------------------------------------------------------
plot(x,y1,'r',x,y0,'b');
xlabel('magnetic field [mT]');
axis tight
legend('12C only','12C+13C');
