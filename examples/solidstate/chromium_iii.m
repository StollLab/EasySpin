% Cr(III) at Q band
%==========================================================================
clear, clf, clc

% parameters and options for pepper
Exp.mwFreq = 34;
Exp.Range = [900 1550];
Exp.nPoints = 4096;
Opt.nKnots = [46 1];
Opt.Verbosity = 1;

% system without magnetic Cr nucleus
Sys.S = 3/2;
Sys.g = 1.990;
Sys.D = [3000 750];     % MHz
Sys.lw = 1;             % mT
[B,spec1] = pepper(Sys,Exp,Opt);

% system with 53Cr nucleus, isotopically pure
Sys.Nucs = '53Cr';
Sys.A = [1 2]*360;
[B,spec2] = pepper(Sys,Exp,Opt);

% system with natural-abundance isotope mixture
Sys.Nucs = 'Cr';
[B,spec3] = pepper(Sys,Exp,Opt);

% display
subplot(3,1,1); plot(B,spec1);
axis tight; title('without 53Cr nucleus');
subplot(3,1,2); plot(B,spec2);
axis tight; title('with 53Cr nucleus');
subplot(3,1,3); plot(B,spec3);
axis tight; title('natural isotope mixture');
xlabel('magnetic field [mT]');
