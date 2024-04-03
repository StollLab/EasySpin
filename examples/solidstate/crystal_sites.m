% single-crystal spectra with separate output of site spectra
%================================================================
clear, clc, clf

% Substitutional nitrogen center (P1) in diamond
P1.g = 2.0024;
P1.Nucs = '14N';
P1.A = [81 114];             % MHz
P1.lwpp = 0.03;              % mT

% Experimental parameters
Exp.mwFreq = 9.5;            % GHz
Exp.CenterSweep = [339 10];  % mT

% P1 and crystal orientation, crystal symmetry
ma = 54.7361;                       % magic angle (deg)
Exp.MolFrame = [45 ma 0]*pi/180;    % P1 molecular frame orientation in crystal
Exp.SampleFrame = [10 20 30]*pi/180;  % crystal orientation in spectrometer
Exp.CrystalSymmetry = 'Fd-3m';      % space group of diamond (#227)

% Simulate spectrum, keeping sites separate
Opt.separate = 'sites';
[B,spc] = pepper(P1,Exp,Opt);

% Plot site spectra separately
stackplot(B,spc);
xlabel('magnetic field (mT)');
