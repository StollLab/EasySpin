function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Gaussian pulse with frequency offset
Params.tp = 0.200; % us
Params.Type = 'gaussian';
Params.tFWHM = 0.064; % us
Params.Amplitude = (pi/Params.tFWHM)/(2*pi);
Params.Frequency = 100; % MHz
Params.TimeStep = 0.0001;

t0 = 0:0.0001:Params.tp;
A = gaussian(t0,Params.tp/2,Params.tFWHM);
A = Params.Amplitude*(A/max(A));
f = cos(2*pi*Params.Frequency*t0)+1i*sin(2*pi*Params.Frequency*t0);
IQ0 = A.*f;

[t,IQ] = pulse(Params);

IQ0 = interp1(t0,IQ0,t,'spline');

err(1) = ~areequal(IQ0,IQ,1e-12,'abs');

% Gaussian pulse with truncation
clear Params
dt = 0.001; % us
t0 = -0.300:dt:0.300; % us
A = gaussian(t0,0,0.100);
A = A/max(A);
ind = find(round(abs(A-0.5)*1e5)/1e5==0);
t0 = t0(ind(1):ind(2))-t0(ind(1));
IQ0 = A(ind(1):ind(2));

Params.tp = t0(end); % us
Params.Type = 'gaussian';
Params.trunc = 0.5;
Params.TimeStep = dt;
Params.Amplitude = 1;

[t,IQ] = pulse(Params);

err(2) = ~areequal(IQ0,IQ,1e-12,'abs');

err = any(err);

data = [];