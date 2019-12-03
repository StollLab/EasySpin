function [err,data] = test(opt,olddata)

% Compare pulse() internal and user-defined sech/tanh
%--------------------------------------------------------------------------

% 1st order sech/tanh with frequency offset
Params.tp = 0.200; % us
Params.Type = 'sech/tanh';
Params.Frequency = [120 0]; % MHz
Params.beta = 10.6;
Params.Qcrit = 5;
Params.TimeStep = 0.0005; % us

[t,IQ] = pulse(Params);


dt = 0.0001; % us
t0 = 0:dt:Params.tp;
% Calculate pulse amplitude from flip angle
BW = (Params.Frequency(2)-Params.Frequency(1))/tanh(Params.beta/2);
Amplitude = sqrt((Params.beta*abs(BW)*Params.Qcrit)/(2*pi*2*Params.tp));
% Amplitude modulation: sech
A = sech((Params.beta/Params.tp)*(t0-Params.tp/2));
% Frequency modulation: tanh
f = ((Params.Frequency(2)-Params.Frequency(1))/(2*tanh(Params.beta/2)))*tanh((Params.beta/Params.tp)*(t0-Params.tp/2));
% Phase modulation
phi = ((Params.Frequency(2)-Params.Frequency(1))/(2*tanh(Params.beta/2)))*(Params.tp/Params.beta)*log(cosh((Params.beta/Params.tp)*(t0-Params.tp/2)));
phi = 2*pi*(phi+mean(Params.Frequency)*t0);
% Pulse
IQ0 = Amplitude*A.*exp(1i*phi);

% Same pulse with user IQ input (and resampling of the time axis)
Params_.tp = 0.200; % us
Params_.Type = 'user-defined';
Params_.I = real(IQ0);
Params_.Q = imag(IQ0);
Params_.TimeStep = Params.TimeStep;

[t_,IQ_] = pulse(Params_);

err = ~areequal(IQ,IQ_,1e-12,'abs');

data = [];