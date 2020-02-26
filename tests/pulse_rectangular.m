function ok = test()

% Compare pulse() output pulse shapes with mathematical expressions

% Rectangular pulse
Params.tp = 0.030; % us
Params.Flip = pi;
Params.TimeStep = 0.001; % us

Amplitude = (Params.Flip/Params.tp)/(2*pi);
t0 = 0:Params.TimeStep:Params.tp;
IQ0(1:numel(t0)) = Amplitude;

[t,IQ] = pulse(Params);

ok = areequal(IQ0,IQ,1e-12,'abs');
