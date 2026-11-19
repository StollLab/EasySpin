% Spin echo of a nitroxide (saffron)
%==========================================================================
% This example demonstrates how to simulate the spin echo of a nitroxide
% powder with shaped pulses.

clear, clc, clf

% Define spin system for nitroxide
Sys.S = 1/2;
Sys.g = [2.00906 2.0062 2.0023];
Sys.Nucs = '14N';
Sys.A = [11.5 11.5 95];  % MHz
Exp.Field = 324.9;  % mT

% Use pepper to simulate the frequency-swept spectrum and to make sure
% pulse excitation bands are set appropriately.
SysCW = Sys;
SysCW.lw = 5;  % MHz
[nu,spc_cw] = pepper(SysCW,Exp);
plot(nu,spc_cw)
xlabel('frequency (GHz)')

%%
% Set up pulses
p90.Type = 'sech/tanh';
p90.beta = 5;
p90.trise = 0.030;  % rise time for smoothed edges, µs
p90.tp = 0.200;  % pulse length, µs
p90.Flip = pi/2;  % flip angle, rad
p90.Frequency = [-120 120];  % excitation band, MHz
p90.Phase = pi;

p180 = p90;
p180.tp = 0.100;  % pulse length, µs
p180.Flip = pi;  % flip angle, rad

% Define pulse sequence
tau = 0.500;  % delay between pulses, µs
Exp.Sequence = {p90, tau, p180, tau+p180.tp};

Exp.mwFreq = 9.1;  % excitation carrier frequency, GHz
Exp.DetFreq = 9.1;  % detection frequency, GHz

% Detect echo transient over time window around end time point in sequence
Exp.DetWindow = [-0.10 0.10]; % start and end time, µs

% Powder average
Opt.GridSize = 20;  % increase number for slower, but more converged simulation

%%
[t,Signal] = saffron(Sys,Exp,Opt);

%%
% Plot echo transient and its FT spectrum

subplot(2,1,1)
plot(t,real(Signal),t,imag(Signal));
xlabel('time relative to the end of the sequence (µs)');
ylabel('signal (arb.u.)');
legend('in phase (I)','out-of-phase (Q)')
title('echo transient')

subplot(2,1,2)
spc_fft = abs(fftshift(fft(Signal)));
nu_fft = fdaxis(t)/1e3 + Exp.DetFreq;  % GHz
plot(nu_fft,spc_fft/max(spc_fft),nu,spc_cw/max(spc_cw));
xlim(Exp.DetFreq+[-1 1]*0.130)
legend('FT of echo','CW EPR');
xlabel('frequency (GHz)')
title('spectrum')
