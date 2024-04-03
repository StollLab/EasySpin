% enhancement of central transition in a S = 3/2 system (spidyan)
%==========================================================================
% this example shows how to enhance the central transition in high-spin
% system using to frequency swept pulses, analogous to 
% Doll, A. et al. Sensitivity enhancement by population transfer in Gd(III)
% spin labels. Phys Chem Chem Phys 17, 7334-7344 (2015).

clear

% System
Sys.S = 3/2;
Sys.ZeemanFreq = 33.500; % GHz
Sys.D = 166;  % MHz

% Pulse
Pulse1.Type = 'quartersin/linear'; % make it a chirp
Pulse1.trise = 0.05; % smooth edges in time domain with a quarter sine, in mus
Pulse1.Qcrit = 10; % use critical adabaticity instead of Par.Flip or Par.Amplitude
Pulse1.tp = 1; % pulse length in µs
Pulse1.Frequency = [500 170]; % frequency band, relative to Exp.mwFreq

Pulse2.Type = 'quartersin/linear'; % make it a chirp
Pulse2.trise = 0.05; % smooth edges in time domain with a quarter sine, in mus
Pulse2.Qcrit = 10; % use critical adabaticity instead of Par.Flip or Par.Amplitude
Pulse2.tp = 1; % pulse length in µs
Pulse2.Frequency = [-500 -170]; % frequency band, relative to Exp.mwFreq

% Sequence
Exp.Sequence = {Pulse1 Pulse2};

Exp.mwFreq = 33.5; % GHz

% The detection operators detect polarization between (1) levels 1 and 2,
% levels 2 and (3) levels 3 and 4
Exp.DetOperator = {'z(1|2)' 'z(2|3)' 'z(3|4)'};

[TimeAxis, Signal] = spidyan(Sys,Exp);

% plotting
figure(1)
clf
plot(TimeAxis,-real(Signal));
xlabel('t ({\mu}s)')
axis tight
ylabel('<S_i>')
legend(Exp.DetOperator)
