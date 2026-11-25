% Selective excitation of two spins with different resonance frequencies (saffron)
%==========================================================================
% shows how to excite two electron spins selectively: the spins have
% different g values, so their resonance frequencies differ by 40 MHz, and
% the narrow-band Gaussian pulses are set to the frequency of one spin
% each, relative to Exp.mwFreq.
%
% With J = 0, the two spins are independent. For an exact selectivity,
% simulate each spin separately (Sys.S = 1/2, with its own g) and apply
% only its own pulses.

clear

% Resonance frequencies of the two spins
nu = [9.48 9.52]; % GHz

Exp.mwFreq = 9.5; % GHz
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.MolFrame = [0 0 0]; % single orientation, molecular frame = lab frame

% System
Sys.S = [1/2 1/2];
Sys.g = gfree*nu(:)/Exp.mwFreq; % g values for the resonance frequencies nu
Sys.J = 0; % MHz

% Selective pulses: Gaussian, excitation width about 18 MHz (FWHM), at the
% offset of spin 1 or 2
offset = (nu - Exp.mwFreq)*1e3; % MHz

P90.Type = 'gaussian';
P90.tFWHM = 0.05; % µs
P90.tp = 0.15; % µs
P90.Flip = pi/2; % rad

P180 = P90;
P180.Flip = pi; % rad

P180_1 = P180; P180_1.Frequency = offset(1);
P90_2 = P90; P90_2.Frequency = offset(2);

tau = 0.5; % µs

% Sequence: pi on spin 1, pi/2 on spin 2, pi on spin 1. Each pulse changes
% the polarization of its own spin only.
Exp.Sequence = {P180_1 tau P90_2 tau P180_1 tau};

Exp.DetWindow = 'sequence'; % detect during the whole sequence
Exp.DetOperator = {'z1' 'z2'};

[TimeAxis, Signal] = saffron(Sys,Exp);

% plotting
figure(1)
clf
plot(TimeAxis,real(Signal));
xlabel('t ({\mu}s)')
axis tight
ylim([-1 1])
ylabel('<S_i>')
legend(Exp.DetOperator)
