function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% 1st order sech/tanh with frequency offset
Exp.tp = 0.200; % us
Exp.PulseShape.Type = 'sech/tanh';
Exp.PulseShape.BW = 120; % MHz
Exp.PulseShape.beta = 10.6;
Exp.PulseShape.SweepDirection = -1;
Exp.Flip = pi;
Exp.TimeStep = 0.0005; % us
Exp.CenterFreq = 60; % MHz

t0 = 0:Exp.TimeStep:Exp.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 5;
BW = Exp.PulseShape.BW/tanh(Exp.PulseShape.beta/2);
Amplitude = sqrt((Exp.PulseShape.beta*BW*Qcrit)/(2*pi*2*Exp.tp));
% Amplitude modulation: sech
A = sech((Exp.PulseShape.beta/Exp.tp)*(t0-Exp.tp/2));
% Frequency modulation: tanh
f = (Exp.PulseShape.BW/(2*tanh(Exp.PulseShape.beta/2)))*tanh((Exp.PulseShape.beta/Exp.tp)*(t0-Exp.tp/2));
% Phase modulation
phi = (Exp.PulseShape.BW/(2*tanh(Exp.PulseShape.beta/2)))*(Exp.tp/Exp.PulseShape.beta)*log(cosh((Exp.PulseShape.beta/Exp.tp)*(t0-Exp.tp/2)));
phi = 2*pi*Exp.PulseShape.SweepDirection*phi;
% Pulse
y0 = Amplitude*A.*exp(1i*(phi+2*pi*Exp.CenterFreq*t0));

Opt.ExciteProfile = 0;
[~,y,~,m] = pulse(Exp,Opt);

suberr(1) = ~areequal(y0,y,1e-12);
suberr(2) = ~areequal(Amplitude*A,m.A,1e-12);
suberr(3) = ~areequal(f,m.nu,1e-12);
suberr(4) = ~areequal(phi,m.phi,1e-12);

err = any(suberr);

data = [];