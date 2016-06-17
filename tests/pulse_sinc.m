function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Sinc pulse
Exp.tp = 0.200; % us
Exp.PulseShape.Type = 'sinc';
Exp.PulseShape.zerocross = 0.034; % us
Exp.TimeStep = 0.001; % us

t0 = 0:Exp.TimeStep:Exp.tp;
A = sin((2*pi*(t0-Exp.tp/2))/Exp.PulseShape.zerocross)./((2*pi*(t0-Exp.tp/2))/Exp.PulseShape.zerocross);
A(t0==Exp.tp/2) = 1;
y0 = A;

[~,y] = pulse(Exp,1);

err = ~areequal(y0,y,1e-12);

data = [];