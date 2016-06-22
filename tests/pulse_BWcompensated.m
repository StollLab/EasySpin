function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------

% Linear chirp and bandwidth-compensated variable rate chirp
Params1.tp = 0.128; % us
Params1.TimeStep = 0.00001; % us
Params1.Type = 'quartersin/linear';
Params1.trise = 0.030; % us
Params1.BW = 500; % MHz
Params1.mwFreq = 9.5; % GHz
Params1.Amplitude = 1; % MHz

Params2 = Params1;
Params2.Type = 'quartersin/BWcompensated';

% Ideal spectrometer magnitude response function
QL = 200; % Q-factor
f0 = 9:0.010:10; % GHz
v1 = abs(1./(1+1i*QL*(f0/Params1.mwFreq-Params1.mwFreq./f0)));

Params2.v1 = v1;
Params2.freqaxis = f0*10^3;

Opt.ExciteProfile = 0;
[t,y1,~,m1] = pulse(Params1,Opt);
[t,y2,~,m2] = pulse(Params2,Opt);

% Calculation
t0 = 0:Params1.TimeStep:Params1.tp;

% Amplitude modulation
A(1:numel(t0)) = 1;
t_part = 0:Params1.TimeStep:Params1.trise;
A(1:numel(t_part)) = sin((t_part)*(pi/(2*Params1.trise)));
A(numel(t0)-numel(t_part)+1:end) = fliplr(A(1:numel(t_part)));
% Frequency modulation for a linear chirp
f = -(Params1.BW/2)+(Params1.BW/Params1.tp)*t0;
% Phase modulation
phi = 2*pi*cumtrapz(t0,f);
phi = phi+abs(min(phi));

% Bandwidth compensation
v1_range = interp1(f0*10^3,v1,f+Params1.mwFreq*10^3);

% Frequency dependence of t and time to frequency mapping
const = trapz(1./v1_range.^2)/t0(end); % const  =  2pi/Qref
t_f = cumtrapz((1/const)*(1./v1_range.^2));
f_adapted = interp1(t_f,f+Params1.mwFreq*10^3,t0,'pchip');
f_adapted = f_adapted-Params1.mwFreq*10^3;
phi_adapted = 2*pi*cumtrapz(t0,f_adapted); % Phase modulation
phi_adapted = phi_adapted+abs(min(phi_adapted));

y0 = A.*exp(1i*phi);
y0_adapted = A.*exp(1i*phi_adapted);

err(1) = ~areequal(y0,y1,1e-12);
err(2) = ~areequal(y0_adapted,y2,1e-12);

err = any(err);

data = [];