function ok = test()

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% 1st order sech/tanh with frequency offset
Params.tp = 0.200; % us
Params.Type = 'sech/tanh';
Params.Frequency = [120 0]; % MHz
Params.beta = 10.6;
Params.Flip = pi;
Params.TimeStep = 0.0005; % us

t0 = 0:Params.TimeStep:Params.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 5;
BW = (Params.Frequency(2)-Params.Frequency(1))/tanh(Params.beta/2);
Amplitude = sqrt((Params.beta*abs(BW)*Qcrit)/(2*pi*2*Params.tp));
% Amplitude modulation: sech
A = sech((Params.beta/Params.tp)*(t0-Params.tp/2));
% Frequency modulation: tanh
f = ((Params.Frequency(2)-Params.Frequency(1))/(2*tanh(Params.beta/2)))*tanh((Params.beta/Params.tp)*(t0-Params.tp/2));
% Phase modulation
phi = ((Params.Frequency(2)-Params.Frequency(1))/(2*tanh(Params.beta/2)))*(Params.tp/Params.beta)*log(cosh((Params.beta/Params.tp)*(t0-Params.tp/2)));
phi = 2*pi*(phi+mean(Params.Frequency)*t0);
% Pulse
IQ0 = Amplitude*A.*exp(1i*phi);

[t,IQ,modulation] = pulse(Params);

ok(1) = areequal(IQ0,IQ,1e-12,'abs');
ok(2) = areequal(Amplitude*A,modulation.A,1e-12,'abs');
ok(3) = areequal(f+mean(Params.Frequency),modulation.freq,1e-12,'abs');
ok(4) = areequal(phi,modulation.phase,1e-12,'abs');
