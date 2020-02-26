function ok = test()

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% HS{6,1}
Params.tp = 0.100; % us
Params.Type = 'sech/uniformQ';
Params.Frequency = [-100 100]; % MHz
Params.beta = 10.4;
Params.n = [6 1];
Params.Amplitude = 30; % MHz
Params.TimeStep = 0.0005; % us

t0 = 0:Params.TimeStep:Params.tp;
npts = numel(t0);
% Amplitude modulation: sech
ti = t0-Params.tp/2;
A(1:round(npts/2)) = sech((Params.beta)*(2^(Params.n(1)-1))*(ti(1:round(npts/2))/Params.tp).^Params.n(1));
A(round(npts/2)+1:npts) = sech((Params.beta)*(2^(Params.n(2)-1))*(ti(round(npts/2)+1:npts)/Params.tp).^Params.n(2));
A = A*Params.Amplitude;
% Frequency modulation
f = cumtrapz(ti,A.^2/trapz(ti,A.^2));
BW = Params.Frequency(2)-Params.Frequency(1);
f = BW*f-BW/2;
% Phase modulation
phi = cumtrapz(ti,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = A.*exp(2i*pi*phi);

[t,IQ,modulation] = pulse(Params);

ok(1) = areequal(IQ0,IQ,1e-11,'abs');
ok(2) = areequal(A,modulation.A,1e-12,'abs');
ok(3) = areequal(f,modulation.freq,1e-12,'abs');
ok(4) = areequal(2*pi*phi,modulation.phase,1e-12,'abs');
