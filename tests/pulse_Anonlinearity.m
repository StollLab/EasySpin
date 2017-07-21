function [err,data] = test(opt,olddata)
% Compare pulse() amplitude compression simulation and amplitude 
% nonlinearity compensation
%--------------------------------------------------------------------------

% Effect of transmitter nonlinearity on pulse shape and amplitude
%--------------------------------------------------------------------------
nu1_max = 15;

% Standard sech amplitude modulation
Params.Type = 'sech/none';
Params.TimeStep = 0.00025; % us
Params.tp = 0.200; % us
Params.beta = 10;
Params.Amplitude = 0.75;

[t,IQ] = pulse(Params);

% 'simulation' effect of amplitude compression
x =  0:0.001:1;
Params.InputAmplitude = x;
Params.OutputAmplitude = nu1_max*(x - 0.30*x.^3); % MHz

Opt.Transmitter = 'simulation';

[t_compressed,IQ_compressed] = pulse(Params,Opt);

[~,ind] = min(abs(x-Params.Amplitude));
suberr(1) = ~areequal(max(abs(IQ_compressed)),Params.OutputAmplitude(ind),1e-11);

% Transmitter nonlinearity compensation test for a sequence of two pulses
%--------------------------------------------------------------------------

clearvars -except suberr

% Two-pulse sequence
%----------------------------------------------
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

% Pulse 1
Params1.Amplitude = amplitudes(1)/max(amplitudes); % relative input amplitude
x =  0:0.001:1;
Params1.InputAmplitude = x;
Params1.OutputAmplitude = nu1_max*(x - 0.30*x.^3); % MHz

Opt.Transmitter = 'compensation';

[t1_compensated,signal1_compensated] = pulse(Params1,Opt);

% Pulse 2
Params2.Amplitude = amplitudes(2)/max(amplitudes); % relative input amplitude
Params2.InputAmplitude = Params1.InputAmplitude;
Params2.OutputAmplitude = Params1.OutputAmplitude; % MHz

[t2_compensated,signal2_compensated] = pulse(Params2,Opt);

% Combined sequence
ParamsTotal.tp = 0.200+0.200+0.100; % us
total_signal = [signal1 zeros(1,0.200/Params1.TimeStep) signal2];
ParamsTotal.I = real(total_signal);
ParamsTotal.Q = imag(total_signal);

ParamsTotal.Amplitude = 1;
ParamsTotal.InputAmplitude = Params1.InputAmplitude;
ParamsTotal.OutputAmplitude = Params1.OutputAmplitude; % MHz

[t_total_compensated,total_signal_compensated] = pulse(ParamsTotal,Opt);

% Compare individually compensated signals with compensated sequence
total_signal = [signal1_compensated zeros(1,0.200/Params1.TimeStep) signal2_compensated];
suberr(2) = ~areequal(total_signal_compensated,total_signal,0.07);

% Effect of transmitter on compensated pulses
%----------------------------------------------
Params1_.tp = 0.200; % us
Params1_.I = real(signal1_compensated);
Params1_.Q = imag(signal1_compensated);

Params1_.Amplitude = max(abs(signal1_compensated))/max(abs(signal2_compensated));
Params1_.InputAmplitude = Params1.InputAmplitude;
Params1_.OutputAmplitude = Params1.OutputAmplitude; % MHz

Opt.Transmitter = 'simulation';

[t1_check,signal1_check] = pulse(Params1_,Opt);

Params2_.tp = 0.100; % us
Params2_.I = real(signal2_compensated);
Params2_.Q = imag(signal2_compensated);

Params2_.Amplitude = max(abs(signal2_compensated))/max(abs(signal2_compensated));
Params2_.InputAmplitude = Params2.InputAmplitude;
Params2_.OutputAmplitude = Params2.OutputAmplitude; % MHz

[t2_check,signal2_check] = pulse(Params2_,Opt);

% Combined sequence (compensated as combined sequence)
ParamsTotal_.tp = 0.200+0.200+0.100; % us
ParamsTotal_.I = real(total_signal_compensated);
ParamsTotal_.Q = imag(total_signal_compensated);

ParamsTotal_.Amplitude = 1;
ParamsTotal_.InputAmplitude = Params2.InputAmplitude;
ParamsTotal_.OutputAmplitude = Params2.OutputAmplitude; % MHz

[t_check,total_signal_check] = pulse(ParamsTotal_,Opt);

% Compare individually compensated signals with compensated sequence
total_signal = [signal1_check zeros(1,0.200/Params1.TimeStep) signal2_check];
suberr(3) = ~areequal(total_signal_check,total_signal,0.07*max(abs(total_signal)));

% Compare compensated and compressed pulse shapes with original pulse shapes
suberr(4) = ~areequal(signal1_check,signal1,0.07*max(abs(signal2)));
suberr(5) = ~areequal(signal2_check,signal2,0.07*max(abs(signal2)));

err = any(suberr);

data = [];