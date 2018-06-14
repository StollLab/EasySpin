clear Exp Sys Opt Pulse
% This script uses transition selective operators to follow the
% polarization build up on the central transition that can be achieved with
% two frequency swept pulses akin to
% Doll, A. et al. Sensitivity enhancement by population transfer in Gd(III)
% spin labels. Phys Chem Chem Phys 17, 7334-7344 (2015).


% System
Sys.S = 3/2;
Sys.ZeemanFreq = 33.500; % GHz
Sys.D = 166;  % MHz

% Pulse
Pulse1.Type = 'quartersin/linear';
Pulse1.trise = 0.05; % us
Pulse1.Qcrit = 10;
Pulse1.tp = 1;
Pulse1.Frequency = [-0.500 -0.17];

Pulse2.Type = 'quartersin/linear';
Pulse2.trise = 0.05; % us
Pulse2.Qcrit = 10;
Pulse2.tp = 1;
Pulse2.Frequency = [0.500 0.17];

% Sequence
Exp.Sequence = {Pulse1 Pulse2};
Exp.Field = 1240; % mT

Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% The detection operators detect polarization between (1) levels 1 and 2,
% levels 2 and (3) levels 3 and 4
Opt.DetOperator = {'z(1|2)' 'z(2|3)' 'z(3|4)'};

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

figure(1)
clf
plot(TimeAxis,-real(Signal));
xlabel('t [{\mu}s]')
axis tight
ylim([-1 3])
ylabel('<S_i>')
legend(Opt.DetOperator)
