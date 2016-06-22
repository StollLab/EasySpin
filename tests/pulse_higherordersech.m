function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% HS8
Params.tp = 0.600; % us
Params.Type = 'sech/uniformQ';
Params.BW = 200; % MHz
Params.beta = 10.6;
Params.n = 8;
Params.Amplitude = 20;
Params.TimeStep = 0.0005; % us

t0 = 0:Params.TimeStep:Params.tp;
% Amplitude modulation: sech
ti = t0-Params.tp/2;
A = sech((Params.beta)*(2^(Params.n-1))*(ti/Params.tp).^Params.n);
A = A*Params.Amplitude;
% Frequency modulation
f = cumtrapz(ti,A.^2/trapz(ti,A.^2));
f = Params.BW*f-Params.BW/2;
% Phase modulation
phi = cumtrapz(ti,f);
phi = phi+abs(min(phi));
% Pulse
y0 = A.*exp(2i*pi*phi);

Opt.ExciteProfile = 0;
[~,y,~,m] = pulse(Params,Opt);

suberr(1) = ~areequal(y0,y,1e-11);
suberr(2) = ~areequal(A,m.A,1e-12);
suberr(3) = ~areequal(f,m.nu,1e-12);
suberr(4) = ~areequal(2*pi*phi,m.phi,1e-12);

err = any(suberr);

data = [];