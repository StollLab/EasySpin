% selective excitation of two spins with identical resonance frequency  (spidyan)
%==========================================================================
% shows how to use spin selective excitation operators. both spins in this
% example have an identical frequency, but by using Opt.ExcOperator it is
% possible to flip them selectively

clear

% System
Sys.S = [1/2 1/2];
Sys.ZeemanFreq = [9.500 9.500]; % GHz

Sys.J = 0; % MHz

P90.Type = 'rectangular';
P90.tp = 0.016; % µs
P90.Flip = pi/2; % rad

P180.Type = 'rectangular';
P180.tp = 0.032; % µs
P180.Flip = pi; % rad

tau = 0.5; % µs

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
