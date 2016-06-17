function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% HS8
Exp.tp = 0.600; % us
Exp.PulseShape.Type = 'sech/uniformQ';
Exp.PulseShape.BW = 200; % MHz
Exp.PulseShape.beta = 10.6;
Exp.PulseShape.n = 8;
Exp.Flip = pi;
Exp.TimeStep = 0.0005; % us

t0 = 0:Exp.TimeStep:Exp.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 8;
Amplitude = sqrt((Exp.PulseShape.beta*Exp.PulseShape.BW*Qcrit)/(2*pi));
% Amplitude modulation: sech
ti = t0-Exp.tp/2;
A = sech((Exp.PulseShape.beta)*(2^(Exp.PulseShape.n-1))*(ti/Exp.tp).^Exp.PulseShape.n);
% Frequency modulation
f = cumtrapz(ti,A.^2/trapz(ti,A.^2));
f = Exp.PulseShape.BW*f-Exp.PulseShape.BW/2;
% Phase modulation
phi = cumtrapz(ti,f);
phi = phi+abs(min(phi));
% Pulse
y0 = Amplitude*A.*exp(2i*pi*phi);

Opt.ExciteProfile = 0;
[~,y,~,m] = pulse(Exp,Opt);

suberr(1) = ~areequal(y0,y,1e-11);
suberr(2) = ~areequal(Amplitude*A,m.A,1e-12);
suberr(3) = ~areequal(f,m.nu,1e-12);
suberr(4) = ~areequal(2*pi*phi,m.phi,1e-12);

err = any(suberr);

data = [];