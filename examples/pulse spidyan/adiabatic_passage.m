clear Exp Sys Opt Pulse

% Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.550; % GHz

% Pulse Definition
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.Amplitude = 20;
Pulse.tp = 0.2;
Pulse.Frequency = [-100 100];
% Pulse.tp=0.01;
% Pulse.Amplitude = 60;

% Experiment/Sequence
Exp.Sequence = {Pulse}; % us
Exp.Field = 1240; % mT
% Exp.TimeStep = 0.00001; % us
 % Frequency range of pulse 

Exp.mwFreq = 33.5000; % GHz

Opt = [];
Opt.SimFreq = 31;

% Options
Exp.DetOperator = {'z1','+1'};
Exp.DetFreq = [0 33.5]; % GHz

% Run simulation
[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting
figure(1)
clf
plot(TimeAxis*1000,real(Signal));
xlabel('t [ns]')
axis tight
ylim([-1 1])
ylabel('<S_i>')
legend(Exp.DetOperator)

% figure(2)
% clf
% plot3(real(Signal(:,2)),imag(Signal(:,2)),real(Signal(:,1)));
% xlabel('<S_x>')
% ylabel('<S_y>')
% zlabel('<S_z>')