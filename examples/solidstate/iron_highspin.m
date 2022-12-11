% High-spin Fe(III), X band CW EPR
%==========================================================================
clear, clc

% Experimental parameters
Exp.mwFreq = 9.5;      % GHz
Exp.Range = [0 700];   % mT - wide sweep
Exp.Temperature = 10;  % K - don't forget to include the temperature!
Exp.nPoints = 5e3;     % lots of points for a wide sweep

% Approach 1: full S=5/2 system with zero-field spliting and isotropic g
Fe_full.S = 5/2;            % 5 unpaired electrons
Fe_full.g = 2;              % isotropic g
Fe_full.D = [10 0.5]*30e3;  % cm^-1 -> MHz
Fe_full.lwpp = 5;           % mT

[B,spc1] = pepper(Fe_full,Exp);

% Approach 2: effective S=1/2 system with anisotropic g (only works for |D|>>kB*T)
Fe_eff.S = 1/2;
Fe_eff.g = [7.13 4.78 1.919];
Fe_eff.lwpp = Fe_full.lwpp;

[B,spc2] = pepper(Fe_eff,Exp);

% Plotting
plot(B,spc1,B,spc2);
grid on
legend('S = 5/2','S_{eff} = 1/2');
legend boxoff
xlabel('magnetic field (mT)');
ylabel('intensity (arb.u.)');
