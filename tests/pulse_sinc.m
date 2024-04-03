function ok = test()

% Compare pulse() output pulse shape with explicit expressions

% Sinc pulse
Params.tp = 0.200;  % µs
Params.Type = 'sinc';
Params.zerocross = 0.034;  % µs
Params.TimeStep = 0.001;  % µs
Params.Amplitude = 1;
[~,IQ] = pulse(Params);

t0 = 0:Params.TimeStep:Params.tp;
dt = t0-Params.tp/2;
x = 2*pi*dt/Params.zerocross;
A = sin(x)./x;
A(dt==0) = 1;
IQref = A;

ok = areequal(IQ,IQref,1e-12,'abs');
