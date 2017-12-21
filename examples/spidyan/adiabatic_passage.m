clear Exp Sys Opt Pulse

% Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

% Pulse Definition
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.Qcrit = 10;

% Experiment/Sequence
Exp.t = 0.2; % us
Exp.Pulses = {Pulse};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100]; % Frequency range of pulse 
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Options
Opt.DetOperator = {'z1','+1'};
Opt.FreqTranslation = [0 -33.5]; % GHz
Opt.FrameShift = 32; % GHz

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
legend(Opt.DetOperator)

figure(2)
clf
plot3(real(Signal(2,:)),imag(Signal(2,:)),real(Signal(1,:)));
xlabel('<S_x>')
ylabel('<S_y>')
zlabel('<S_z>')