function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Linear chirp with frequency offset
Exp.tp = 0.064; % us
Exp.PulseShape.Type = 'rectangular/linear';
Exp.PulseShape.BW = 120; % MHz
Exp.Flip = pi/2;
Exp.TimeStep = 0.0001; % us
Exp.CenterFreq = 120; % MHz

t0 = 0:Exp.TimeStep:Exp.tp; % us
% Pulse amplitude (Qcrit = 2ln2/pi for a pi/2 pulse)
Amplitude = sqrt((4*log(2)*Exp.PulseShape.BW)/(Exp.tp))/(2*pi);
% Frequency modulation: linear chirp
f = -(Exp.PulseShape.BW/2)+(Exp.PulseShape.BW/Exp.tp)*t0;
% Phase modulation
phi = cumtrapz(t0,f);
phi = phi+abs(min(phi));
% Pulse
y0 = Amplitude*exp(2i*pi*(phi+Exp.CenterFreq*t0));

[~,y] = pulse(Exp);

err(1) = ~areequal(y0,y,1e-12);

% Quartersin weighted chirp
clearvars -except err
Exp.tp = 0.128; % us
Exp.PulseShape.Type = 'quartersin/linear';
Exp.PulseShape.trise = 0.010; % us
Exp.PulseShape.BW = 100; % MHz
Exp.Amplitude = 20; % MHz
Opt.OverSampleFactor = 8;

dt = 1/(2*Opt.OverSampleFactor*Exp.PulseShape.BW/2);
dt = Exp.tp/(round(Exp.tp/dt));
t0 = 0:dt:Exp.tp; % us
% Amplitude modulation: quarter sine weighting
t_part = 0:dt:Exp.PulseShape.trise;
A(1:numel(t0)) = 1;
A(1:numel(t_part)) = sin((pi*t_part)/(2*Exp.PulseShape.trise));
A(end-numel(t_part)+1:end) = fliplr(A(1:numel(t_part)));
% Frequency modulation: linear chirp
f = -(Exp.PulseShape.BW/2)+(Exp.PulseShape.BW/Exp.tp)*t0;
% Phase modulation
phi = cumtrapz(t0,f);
phi = phi+abs(min(phi));
% Pulse
y0 = Exp.Amplitude*A.*exp(2i*pi*phi);

Opt.ExciteProfile = 0;
[~,y,~,m] = pulse(Exp,Opt);

suberr(1) = ~areequal(y0,y,1e-12);
suberr(2) = ~areequal(Exp.Amplitude*A,m.A,1e-12);
suberr(3) = ~areequal(f,m.nu,1e-12);

err(2) = any(suberr);

err = any(err);

data = [];