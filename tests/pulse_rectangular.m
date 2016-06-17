function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Rectangular pulse
Exp.tp = 0.030; % us
Exp.Flip = pi;
Exp.TimeStep = 0.001; % us

Amplitude = (Exp.Flip/Exp.tp)/(2*pi);
t0 = 0:Exp.TimeStep:Exp.tp;
y0(1:numel(t0)) = Amplitude;

[~,y] = pulse(Exp);

err = ~areequal(y0,y,1e-12);

data = [];
