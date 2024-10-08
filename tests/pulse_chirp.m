function ok = test()

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Linear chirp with frequency offset
Params.tp = 0.064;  % µs
Params.Type = 'rectangular/linear';
Params.Frequency = [-60 60] + 120;  % MHz
Params.Flip = pi/2;
Params.TimeStep = 0.0001;  % µs

t0 = 0:Params.TimeStep:Params.tp;  % µs
% Pulse amplitude (Qcrit = 2ln2/pi for a pi/2 pulse)
BW = Params.Frequency(2)-Params.Frequency(1);
Amplitude = sqrt((4*log(2)*BW/Params.tp))/(2*pi);
% Frequency modulation: linear chirp
f = -(BW/2)+(BW/Params.tp)*t0;
% Phase modulation
phi = cumtrapz(t0,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = Amplitude*exp(2i*pi*(phi+mean(Params.Frequency)*t0));

[~,IQ] = pulse(Params);

ok(1) = areequal(IQ0,IQ,1e-12,'abs');

% Quartersin weighted chirp
clear Params
Params.tp = 0.128;  % µs
Params.Type = 'quartersin/linear';
Params.trise = 0.010;  % µs
Params.Frequency = [-50 50];  % MHz
Params.Amplitude = 20;  % MHz
Opt.OverSampleFactor = 10;

BW = Params.Frequency(2)-Params.Frequency(1);
dt = 1/(2*Opt.OverSampleFactor*BW/2);
dt = Params.tp/(round(Params.tp/dt));
t0 = 0:dt:Params.tp; % µs
% Amplitude modulation: quarter sine weighting
t_part = 0:dt:Params.trise;
A(1:numel(t0)) = 1;
A(1:numel(t_part)) = sin((pi*t_part)/(2*Params.trise));
A(end-numel(t_part)+1:end) = fliplr(A(1:numel(t_part)));
% Frequency modulation: linear chirp
f = -(BW/2)+(BW/Params.tp)*t0;
% Phase modulation
phi = cumtrapz(t0,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = Params.Amplitude*A.*exp(2i*pi*phi);

[~,IQ,modulation] = pulse(Params);

ok(2) = areequal(IQ0,IQ,1e-12,'abs');
ok(3) = areequal(Params.Amplitude*A,modulation.A,1e-12,'abs');
ok(4) = areequal(f,modulation.freq,1e-12,'abs');


% Halfsin weighted chirp
clear Params
Params.tp = 0.128;  % µs
Params.Type = 'halfsin/linear';
Params.Frequency = [-50 50];  % MHz
Params.Amplitude = 20;  % MHz
Opt.OverSampleFactor = 10;

BW = Params.Frequency(2)-Params.Frequency(1);
dt = 1/(2*Opt.OverSampleFactor*BW/2);
dt = Params.tp/(round(Params.tp/dt));
t0 = 0:dt:Params.tp; % µs
% Amplitude modulation: half sine weighting
A = sin(pi*t0/Params.tp);
% Frequency modulation: linear chirp
f = -(BW/2)+(BW/Params.tp)*t0;
% Phase modulation
phi = cumtrapz(t0,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = Params.Amplitude*A.*exp(2i*pi*phi);

[~,IQ,modulation] = pulse(Params);

ok(5) = areequal(IQ0,IQ,1e-12,'abs');
ok(6) = areequal(Params.Amplitude*A,modulation.A,1e-12,'abs');
ok(7) = areequal(f,modulation.freq,1e-12,'abs');
