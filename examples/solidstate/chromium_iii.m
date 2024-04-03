% Cr(III) at Q band
%==========================================================================
clear, clf, clc

% Experimental parameters
Exp.mwFreq = 34;
Exp.Range = [900 1550];
Exp.nPoints = 4096;

% Simulation options
Opt.GridSize = [46 1];
Opt.Verbosity = 1;

% System without magnetic Cr nucleus
Sys.S = 3/2;
Sys.g = 1.990;
Sys.D = [3000 750];     % MHz
Sys.lw = 1;             % mT
[B,spec0] = pepper(Sys,Exp,Opt);

% System with 53Cr nucleus, isotopically pure
Sys.Nucs = '53Cr';
Sys.A = [1 2]*360;  % MHz
[B,spec_53] = pepper(Sys,Exp,Opt);

% System with natural-abundance isotope mixture
Sys.Nucs = 'Cr';
[B,spec_na] = pepper(Sys,Exp,Opt);

% Plotting
subplot(3,1,1);
plot(B,spec0);
axis tight;
title('without 53Cr nucleus');
subplot(3,1,2);
plot(B,spec_53);
axis tight;
title('with 53Cr nucleus');
subplot(3,1,3);
plot(B,spec_na);
axis tight;
title('natural isotope mixture');
xlabel('magnetic field (mT)');
