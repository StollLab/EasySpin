% Phenalenyl radical anion, including 13C satellite lines
%==========================================================================
% Here, we manually simulate the spectra of each 12C and 13C
% isotopologue of the phenalenyl radical and sum them to get 
% the total spectrum corresponding to natural-abundance carbon.
% (We neglect species with more than one 13C.)
%
% EasySpin can do this automatically, see phenalenyl.m
%
% Hyperfine values taken from Gerson book, Table 8.4, page 225

clear, clf

% Spin system
%--------------------------------------------------------------------------
A_H = unitconvert([-0.629 +0.181],'mT->MHz');
A_C = unitconvert([+0.966 -0.784 -0.784],'mT->MHz');
Sys.g = 2;
Sys.Nucs = '1H,1H';
Sys.n = [6 3];
Sys.A = A_H;
Sys.lwpp = [0,0.01];  % Lorentzian line shape

% Experimental parameters
%--------------------------------------------------------------------------
Exp.mwFreq = 9.5;
Exp.Range = [336.5 342.5];
Exp.nPoints = 10000;

% Simulations
%--------------------------------------------------------------------------

% Spectrum of species with all 12C 
[B,spc0] = garlic(Sys,Exp);

% Spectra with species containing one 13C
Sys.Nucs = '1H,1H,13C';
Sys.n = [6 3 1];
% 13C in one of positions 1,3,4,6,7,9
Sys.A = [A_H A_C(1)];
[B,spc1] = garlic(Sys,Exp);
% 13C in one of positions 2,5,8
Sys.A = [A_H A_C(2)];
[B,spc2] = garlic(Sys,Exp);
% 13C in one of positions 3a,6a,9a
Sys.A = [A_H A_C(3)];
[B,spc3] = garlic(Sys,Exp);

abund = nucabund('12C,13C');
aC12 = abund(1)^12;
aC13 = abund(1)^11*abund(2)^1;
y_total = aC12*spc0 + aC13*(6*spc1 + 3*spc2 + 3*spc3);

% Species with more than one 13C are neglected.


% Graphical rendering
%--------------------------------------------------------------------------
plot(B,y_total);
xlabel('magnetic field (mT)');
axis tight
