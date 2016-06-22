function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Rectangular pulse
Params.tp = 0.030; % us
Params.Flip = pi;
Params.TimeStep = 0.001; % us

Amplitude = (Params.Flip/Params.tp)/(2*pi);
t0 = 0:Params.TimeStep:Params.tp;
y0(1:numel(t0)) = Amplitude;

[~,y] = pulse(Params);

err = ~areequal(y0,y,1e-12);

data = [];
