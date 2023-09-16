% single-crystal spectra with separate output of site spectra
%================================================================
clear, clc

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

Opt.separate = 'sites';

[B,spc,info] = pepper(P1,Exp,Opt);
plot(B,spc);
