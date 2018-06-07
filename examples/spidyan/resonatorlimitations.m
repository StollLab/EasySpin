% Simulates three Electron spins with different ZeemanFrequencies during an
% adiabatic pulse with (see comments in for loop)
% - No resonator
% - Resonator
% - And compensation of the resonator

clear Sys Opt Exp Pulse

% Zeeman frequencies of the spins
ZFreqs = [33.4 33.5 33.6];

Sys.S = [1/2];

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.1; % us
Pulse.Qcrit = 5; % If a critical adiabaticity is provided for the Pulse,
% Exp.Flip does not need to be defined or will be ignored

Exp.t = 0.5; % us
Exp.Pulses = {Pulse};
Exp.Frequency = [-0.15 0.15]; % frequency range of sweep, GHz

Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Uncomment the following two lines to see the effect of the resonator
Exp.Resonator.nu0 = 33.5;
Exp.Resonator.QL = 300;

% Uncomment the following line to compensate for the bandwidth limitation
Exp.Resonator.Mode = 'compensate';

Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;


% Set up for plotting
figure(1)
clf
hold on

Signal = cell(1,length(ZFreqs));

for i = 1 : length(ZFreqs)
  
  Sys.ZeemanFreq = [ZFreqs(i)];
  
  [TimeAxis, Signal{i}] = spidyan(Sys,Exp,Opt);
  
  plot(TimeAxis*1000,real(Signal{i}));
  
end

xlabel('t [ns]')
axis tight
ylim([-1 1])
legend(string(ZFreqs))