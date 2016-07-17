% comparison of field sweep and frequency sweep cw EPR
%==========================================================================
clear, clf

% Spin system: nitroxide radical
%------------------------------------------------------------
Nx.g = [2.008 2.006 2.002];
Nx.Nucs = '14N';
Nx.A = [20 20 90]; % MHz


% Frequency-sweep experiment
%------------------------------------------------------------
% For frequency sweeps, the convolutional width should be given in MHz.
Nx.lwpp = 5; % MHz
Exp1.Field = 3300; % static magnetic field, mT
Exp1.mwRange = [92 93]; % frequency sweep range, GHz

% Simulate frequency-sweep spectrum
[nu,specnu] = pepper(Nx,Exp1);

% Field-sweep experiment
%------------------------------------------------------------
% For field sweeps, the convolutional width should be given in mT.
Nx.lwpp = 0.2; % mT
Exp2.mwFreq = 92.6; % constant microwave frequency, GHz
Exp2.Range = [3285 3315]; % field sweep range, mT
Exp2.Harmonic = 0;

% Simulate field-sweep spectrum
[B,specB] = pepper(Nx,Exp2);

% Plot the results
%------------------------------------------------------------
subplot(2,1,1);
plot(nu,specnu);
title('frequency sweep');
xlabel('frequency (GHz)');

subplot(2,1,2);
plot(B/1e3,specB);
title('field sweep');
xlabel('magnetic field (T)');
