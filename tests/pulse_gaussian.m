function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Gaussian pulse with frequency offset
Exp.tp = 0.200; % us
Exp.PulseShape.Type = 'gaussian';
Exp.PulseShape.tFWHM = 0.064; % us
Exp.Amplitude = (pi/Exp.PulseShape.tFWHM)/(2*pi);
Exp.CenterFreq = 100; % MHz

t0 = 0:0.0001:Exp.tp;
A = gaussian(t0,Exp.tp/2,Exp.PulseShape.tFWHM);
A = Exp.Amplitude*(A/max(A));
f = cos(2*pi*Exp.CenterFreq*t0)+1i*sin(2*pi*Exp.CenterFreq*t0);
y0 = A.*f;

[t,y] = pulse(Exp);

y0 = interp1(t0,y0,t,'spline');

err(1) = ~areequal(y0,y,1e-12);

% Gaussian pulse with truncation
clearvars -except err
dt = 0.001; % us
t0 = -0.300:dt:0.300; % us
A = gaussian(t0,0,0.100);
A = A/max(A);
ind = find(round(abs(A-0.5)*1e5)/1e5==0);
t0 = t0(ind(1):ind(2))-t0(ind(1));
y0 = A(ind(1):ind(2));

Exp.tp = t0(end); % us
Exp.PulseShape.Type = 'gaussian';
Exp.PulseShape.trunc = 0.5;
Exp.TimeStep = dt;

[~,y] = pulse(Exp);

err(2) = ~areequal(y0,y,1e-12);

err = any(err);

data = [];