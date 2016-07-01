function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% 1st order sech/tanh with frequency offset
Params.tp = 0.200; % us
Params.Type = 'sech/tanh';
Params.BW = 120; % MHz
Params.beta = 10.6;
Params.SweepDirection = -1;
Params.Flip = pi;
Params.TimeStep = 0.0005; % us
Params.CenterFreq = 60; % MHz

t0 = 0:Params.TimeStep:Params.tp;
% Calculate pulse amplitude from flip angle
Qcrit = 5;
BW = Params.BW/tanh(Params.beta/2);
Amplitude = sqrt((Params.beta*BW*Qcrit)/(2*pi*2*Params.tp));
% Amplitude modulation: sech
A = sech((Params.beta/Params.tp)*(t0-Params.tp/2));
% Frequency modulation: tanh
f = (Params.BW/(2*tanh(Params.beta/2)))*tanh((Params.beta/Params.tp)*(t0-Params.tp/2));
% Phase modulation
phi = (Params.BW/(2*tanh(Params.beta/2)))*(Params.tp/Params.beta)*log(cosh((Params.beta/Params.tp)*(t0-Params.tp/2)));
phi = 2*pi*Params.SweepDirection*phi;
% Pulse
IQ0 = Amplitude*A.*exp(1i*(phi+2*pi*Params.CenterFreq*t0));

Opt.ExciteProfile = 0;
[t,IQ,exprofile,modulation] = pulse(Params,Opt);

suberr(1) = ~areequal(IQ0,IQ,1e-12);
suberr(2) = ~areequal(Amplitude*A,modulation.A,1e-12);
suberr(3) = ~areequal(f+Params.CenterFreq,modulation.nu,1e-12);
suberr(4) = ~areequal(phi,modulation.phi,1e-12);

err = any(suberr);

data = [];