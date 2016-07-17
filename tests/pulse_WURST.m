function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% HS8
Params.tp = 0.500; % us
Params.Type = 'WURST/linear';
Params.Frequency = [-250 250]+100; % MHz
Params.nwurst = 15;
Params.Amplitude = 15; % MHz

Opt.ExciteProfile = false;
[t,IQ,exprofile,modulation] = pulse(Params,Opt);

% Amplitude modulation: WURST
A = 1-abs((sin((pi*(t-Params.tp/2))/Params.tp)).^Params.nwurst);
% Frequency modulation
BW = Params.Frequency(2)-Params.Frequency(1);
f = -(BW/2)+(BW/Params.tp)*t;
% Phase modulation
phi = cumtrapz(t-Params.tp/2,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = Params.Amplitude*A.*exp(2i*pi*(phi+mean(Params.Frequency)*t));

suberr(1) = ~areequal(IQ0,IQ,1e-11);
suberr(2) = ~areequal(Params.Amplitude*A,modulation.A,1e-12);
suberr(3) = ~areequal(f+mean(Params.Frequency),modulation.freq,1e-12);
suberr(4) = ~areequal(2*pi*phi,modulation.phase,1e-12);

err = any(suberr);

data = [];