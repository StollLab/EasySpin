function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% HS8
Params.tp = 0.500; % us
Params.Type = 'WURST/linear';
Params.BW = 500; % MHz
Params.nwurst = 15;
Params.Amplitude = 15; % MHz
Params.CenterFreq = 100; % MHz

Opt.ExciteProfile = 0;
[t,IQ,exprofile,modulation] = pulse(Params,Opt);

% Amplitude modulation: WURST
A = 1-abs((sin((pi*(t-Params.tp/2))/Params.tp)).^Params.nwurst);
% Frequency modulation
f = -(Params.BW/2)+(Params.BW/Params.tp)*t;
% Phase modulation
phi = cumtrapz(t-Params.tp/2,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = Params.Amplitude*A.*exp(2i*pi*(phi+Params.CenterFreq*t));

suberr(1) = ~areequal(IQ0,IQ,1e-11);
suberr(2) = ~areequal(Params.Amplitude*A,modulation.A,1e-12);
suberr(3) = ~areequal(f+Params.CenterFreq,modulation.nu,1e-12);
suberr(4) = ~areequal(2*pi*phi,modulation.phi,1e-12);

err = any(suberr);

data = [];