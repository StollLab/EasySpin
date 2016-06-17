function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% HS8
Exp.tp = 0.500; % us
Exp.PulseShape.Type = 'WURST/linear';
Exp.PulseShape.BW = 500; % MHz
Exp.PulseShape.n = 15;
Exp.Amplitude = 15; % MHz
Exp.CenterFreq = 100; % MHz

Opt.ExciteProfile = 0;
[t,y,~,m] = pulse(Exp,Opt);

% Amplitude modulation: WURST
A = 1-abs((sin((pi*(t-Exp.tp/2))/Exp.tp)).^Exp.PulseShape.n);
% Frequency modulation
f = -(Exp.PulseShape.BW/2)+(Exp.PulseShape.BW/Exp.tp)*t;
% Phase modulation
phi = cumtrapz(t-Exp.tp/2,f);
phi = phi+abs(min(phi));
% Pulse
y0 = Exp.Amplitude*A.*exp(2i*pi*(phi+Exp.CenterFreq*t));

suberr(1) = ~areequal(y0,y,1e-11);
suberr(2) = ~areequal(Exp.Amplitude*A,m.A,1e-12);
suberr(3) = ~areequal(f,m.nu,1e-12);
suberr(4) = ~areequal(2*pi*phi,m.phi,1e-12);

err = any(suberr);

data = [];