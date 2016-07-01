function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Linear chirp with frequency offset
Params.tp = 0.064; % us
Params.Type = 'rectangular/linear';
Params.BW = 120; % MHz
Params.Flip = pi/2;
Params.TimeStep = 0.0001; % us
Params.CenterFreq = 120; % MHz

t0 = 0:Params.TimeStep:Params.tp; % us
% Pulse amplitude (Qcrit = 2ln2/pi for a pi/2 pulse)
Amplitude = sqrt((4*log(2)*Params.BW/Params.tp))/(2*pi);
% Frequency modulation: linear chirp
f = -(Params.BW/2)+(Params.BW/Params.tp)*t0;
% Phase modulation
phi = cumtrapz(t0,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = Amplitude*exp(2i*pi*(phi+Params.CenterFreq*t0));

[t,IQ] = pulse(Params);

err(1) = ~areequal(IQ0,IQ,1e-12);

% Quartersin weighted chirp
clearvars -except err
Params.tp = 0.128; % us
Params.Type = 'quartersin/linear';
Params.trise = 0.010; % us
Params.BW = 100; % MHz
Params.Amplitude = 20; % MHz
Opt.OverSampleFactor = 8;

dt = 1/(2*Opt.OverSampleFactor*Params.BW/2);
dt = Params.tp/(round(Params.tp/dt));
t0 = 0:dt:Params.tp; % us
% Amplitude modulation: quarter sine weighting
t_part = 0:dt:Params.trise;
A(1:numel(t0)) = 1;
A(1:numel(t_part)) = sin((pi*t_part)/(2*Params.trise));
A(end-numel(t_part)+1:end) = fliplr(A(1:numel(t_part)));
% Frequency modulation: linear chirp
f = -(Params.BW/2)+(Params.BW/Params.tp)*t0;
% Phase modulation
phi = cumtrapz(t0,f);
phi = phi+abs(min(phi));
% Pulse
IQ0 = Params.Amplitude*A.*exp(2i*pi*phi);

Opt.ExciteProfile = 0;
[t,IQ,exprofile,modulation] = pulse(Params,Opt);

suberr(1) = ~areequal(IQ0,IQ,1e-12);
suberr(2) = ~areequal(Params.Amplitude*A,modulation.A,1e-12);
suberr(3) = ~areequal(f,modulation.nu,1e-12);

err(2) = any(suberr);

err = any(err);

data = [];