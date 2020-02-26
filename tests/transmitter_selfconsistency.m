function ok = test()

% Compare transmitter() amplitude compression simulation and amplitude 
% nonlinearity compensation (self-consistency check)
%--------------------------------------------------------------------------

% Transmitter nonlinearity compensation test for a sequence of two pulses
%--------------------------------------------------------------------------
amplitudes = [13 60]; % MHz

% Pulse 1
Params1.Type = 'sech/tanh';
Params1.TimeStep = 0.00025; % us
Params1.Frequency = [-50 50]; % MHz
Params1.tp = 0.200; % us
Params1.beta = 10;
Params1.Amplitude = 13; % MHz

[t1,signal1] = pulse(Params1);

% Pulse 2
Params2 = Params1;
Params2.tp = 0.100; % us
Params2.Amplitude = 60; % MHz
[t2,signal2] = pulse(Params2);


% Amplitude compensation
%----------------------------------------------
nu1_max = 85; % MHz
nu1_out = 50; % MHz

Ain =  0:0.001:1; % relative scale
Aout =  nu1_max*(Ain - 0.30*Ain.^3); % MHz

% Pulse 1
signal1_compensated = transmitter(signal1*nu1_out/max(amplitudes),Ain,Aout,'compensate');

% Pulse 2
signal2_compensated = transmitter(signal2*nu1_out/max(amplitudes),Ain,Aout,'compensate');

% Combined sequence
total_signal = [signal1 zeros(1,0.200/Params1.TimeStep) signal2];

total_signal_compensated = transmitter(total_signal*nu1_out/max(amplitudes),Ain,Aout,'compensate');

% Compare individually compensated signals with compensated sequence
total_signal_combined = [signal1_compensated zeros(1,0.200/Params1.TimeStep) signal2_compensated];
ok(1) = areequal(total_signal_compensated,total_signal_combined,1e-6,'abs');


% Effect of transmitter on compensated pulses
%----------------------------------------------
% Pulse 1
signal1_check = transmitter(signal1_compensated,Ain,Aout,'simulate');

% Pulse 2
signal2_check = transmitter(signal2_compensated,Ain,Aout,'simulate');

% Combined sequence (compensated as combined sequence)
total_signal_check = transmitter(total_signal_compensated,Ain,Aout,'simulate');

% Compare individually compensated signals with compensated sequence
total_signal_combined = [signal1_check zeros(1,0.200/Params1.TimeStep) signal2_check];
ok(2) = areequal(total_signal_check,total_signal_combined,1e-6,'abs');

% Compare compensated and compressed pulse shapes with original pulse shapes
scalefactor = nu1_out/amplitudes(2);
ok(3) = areequal(signal1_check,scalefactor*signal1,1e-6,'abs');
ok(4) = areequal(signal2_check,scalefactor*signal2,1e-6,'abs');
