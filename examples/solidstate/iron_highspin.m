% High-spin Fe(III), X band CW EPR
%==========================================================================
clear, clc

% Experimental parameters
Exp.mwFreq = 9.5; % GHz
Exp.Range = [0 800]; % mT - wide sweep
Exp.Temperature = 10; % K- don't forget to include the temperature!
Exp.nPoints = 5e3; % lots of points for a wide sweep

% Approach 1: S=5/2 system with zero-field spliting and isotropic g
Fe1.S = 5/2;
Fe1.g = 2;
cm2MHz = 29e3; % conversion from cm^-1 to MHz
Fe1.D = [1 0.05]*10*cm2MHz; % MHz
Fe1.lwpp = 5; % mT
[B,spc1] = pepper(Fe1,Exp);

% Approach 2: effective S=1/2 system with anisotropic g
% (This only works for |D|>>kB*T)
Fe2.S = 1/2;
Fe2.g = [7.13 4.78 1.919];
Fe2.lwpp = Fe1.lwpp;
[B,spc2] = pepper(Fe2,Exp);

% Plotting
plot(B,spc1,B,spc2);
grid on
legend('S=5/2','S_{eff}=1/2');
xlabel('magnetic field (mT)');
ylabel('intensity (arb.u.)');
