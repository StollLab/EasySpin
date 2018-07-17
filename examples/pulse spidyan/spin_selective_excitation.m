% selective excitation of two spins with identical resonance frequency  (spidyan)
%==========================================================================
% this example shows how to enhance the central transition in high-spin
% system using to frequency swept pulses, analogous to 
% Doll, A. et al. Sensitivity enhancement by population transfer in Gd(III)
% spin labels. Phys Chem Chem Phys 17, 7334-7344 (2015).

clear

% System
Sys.S = [1/2 1/2];
Sys.ZeemanFreq = [9.500 9.500]; % GHz

Sys.J = 0; % MHz

P90.tp = 0.016; % mus
P90.Flip = pi/2; % rad

P180.tp = 0.032; % mus
P180.Flip = pi; % rad

tau = 0.5; % mus

% Sequence
Exp.Sequence = {P90 tau P180 tau P180 tau};

Exp.mwFreq = 9.5; % GHz

Exp.DetOperator = {'z1' 'z2'};

Opt.ExcOperator = {'x1' [] 'x2'}; % spin selective pulses

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% plotting
figure(1)
clf
plot(TimeAxis,real(Signal));
xlabel('t ({\mu}s)')
axis tight
ylim([-1 1])
ylabel('<S_i>')
legend(Exp.DetOperator)
