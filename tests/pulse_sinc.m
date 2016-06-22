function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Sinc pulse
Params.tp = 0.200; % us
Params.Type = 'sinc';
Params.zerocross = 0.034; % us
Params.TimeStep = 0.001; % us
Params.Amplitude = 1;

t0 = 0:Params.TimeStep:Params.tp;
A = sin((2*pi*(t0-Params.tp/2))/Params.zerocross)./((2*pi*(t0-Params.tp/2))/Params.zerocross);
A(t0==Params.tp/2) = 1;
y0 = A;

[~,y] = pulse(Params);

err = ~areequal(y0,y,1e-12);

data = [];