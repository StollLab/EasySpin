% resonator limitations and compensation (spidyan)
%==========================================================================
% this script shows the effect of a resonator on in the inversion of three
% spins. One spin has a resonance frequency that is the same as the
% centerfrequency of the resonator, one is +100 MHz from the resonator
% frequency and the third is -100 MHz from the resonator center frequency
% You can have a look at these three options:
% - No resonator
% - Resonator
% - And compensation of the resonator
%
% Compensation of the resonator profile was described in 
% Doll, A., Pribitzer, S., Tschaggelar, R., Jeschke, G.,
% J. Magn. Reson. 230, 27-39 (2013), DOI: 10.1016/j.jmr.2013.01.002

clear

% Zeeman frequencies of the spins
ZFreqs = [9.4 9.5 9.6]; % GHz

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.1; % µs
Pulse.Qcrit = 5; % If a critical adiabaticity is provided for the Pulse,
% Pulse.Flip does not need to be defined or will be ignored
Pulse.tp = 0.5; % µs
Pulse.Frequency = [-200 200]; % MHz

Exp.Sequence = {Pulse};
Exp.mwFreq = 9.5; % GHz
Exp.DetOperator = {'z1'};

% Uncomment the following two lines to see the effect of the resonator
Exp.ResonatorFrequency = 9.5; % resonator frequency in GHz
Exp.ResonatorQL = 150; % loaded Q factor of resonator

% Uncomment the following line to compensate for the bandwidth limitation,
% also note how the pulses in the plot get longer than 500 ns
% Exp.ResonatorMode = 'compensate';

% Set up for plotting
figure(1)
clf
hold on

Signal = cell(1,length(ZFreqs));

for i = 1 : length(ZFreqs)
  
  Sys_.ZeemanFreq = [ZFreqs(i)];
  [TimeAxis, Signal{i}] = spidyan(Sys_,Exp);
  
  plot(TimeAxis*1000,real(Signal{i}));
  
end

xlabel('t (ns)')
axis tight
ylim([-1 1])
legend(string(ZFreqs))