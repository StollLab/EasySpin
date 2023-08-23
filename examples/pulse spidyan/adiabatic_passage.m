% spin trajectory during adiabatic passage (spidyan)
%==========================================================================
% With this script the trajectory of spin during adiabatic passage is
% plotted in Fig. 1 in two and in Fig 2 in three dimensions

clear

% Spin System
Sys.ZeemanFreq = 9.500; % resonance frequency of the spin in GHz

% Pulse Definition
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.Qcrit = 5; % critical adiabaticity
Pulse.tp = 0.2; % µs
Pulse.Frequency = [-100 100]; % frequency band, MHz

% Experiment/Sequence
Exp.Sequence = {Pulse}; % µs
Exp.mwFreq = 9.500; % GHz

% Options
Exp.DetOperator = {'z1','+1'}; % detection operators
Exp.DetFreq = [0 9.5]; % downconversion frequency in GHz

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp);

% Plotting
figure(1)
clf
plot(TimeAxis*1000,real(Signal));
xlabel('t (ns)')
axis tight
ylim([-1 1])
ylabel('<S_i>')
legend(Exp.DetOperator)

figure(2)
clf
plot3(real(Signal(:,2)),imag(Signal(:,2)),real(Signal(:,1)));
xlabel('<S_x>')
ylabel('<S_y>')
zlabel('<S_z>')
grid on