function [err,data] = test(opt,olddata)

% Compare pulse() output pulse shapes with mathematical expressions
%--------------------------------------------------------------------------
%%
% Linear chirp and bandwidth-compensated variable rate chirp
%--------------------------------------------------------------------------
Params.tp = 0.128; % us
Params.TimeStep = 0.00001; % us
Params.Type = 'quartersin/linear';
Params.trise = 0.030; % us
Params.Frequency = [-150 150]; % MHz
Params.Amplitude = 1; % MHz

Params.mwFreq = 9.5; % GHz

% Ideal spectrometer magnitude response function
Params.ResonatorFrequency = Params.mwFreq;
Params.ResonatorQL = 50; % Q-factor

[t,IQ1] = pulse(Params);
Opt.Resonator = 'BWcompensation';
[t,IQ2] = pulse(Params,Opt);

% Consistency check in the frequency domain
f = fdaxis(t);
IQ1_ft = abs(fftshift(fft(IQ1)));
IQ1_ft = IQ1_ft/max(IQ1_ft);
IQ2_ft = abs(fftshift(fft(IQ2)));
IQ2_ft = IQ2_ft/max(IQ2_ft);

f0 = 9:0.00001:10; % GHz
H = 1./(1+1i*Params.ResonatorQL*(f0/Params.ResonatorFrequency-Params.ResonatorFrequency./f0));
v1 = abs(H);
profile = interp1((f0-Params.mwFreq)*1e3,v1,f);
profile = profile/max(profile);
check = profile.*IQ2_ft;
check = check/max(check);

err(1) = ~areequal(std(IQ1_ft(abs(f)<0.6*Params.Frequency(2))),std(check(abs(f)<0.6*Params.Frequency(2))),1e-2);

% Calculation
t0 = 0:Params.TimeStep:Params.tp;

% Amplitude modulation
A(1:numel(t0)) = 1;
t_part = 0:Params.TimeStep:Params.trise;
A(1:numel(t_part)) = sin((t_part)*(pi/(2*Params.trise)));
A(numel(t0)-numel(t_part)+1:end) = fliplr(A(1:numel(t_part)));
% Frequency modulation for a linear chirp
BW = Params.Frequency(2)-Params.Frequency(1);
f = -(BW/2)+(BW/Params.tp)*t0;
% Phase modulation
phi = 2*pi*cumtrapz(t0,f);
phi = phi+abs(min(phi));

% Bandwidth compensation
v1_range = interp1(f0*10^3,v1,f+Params.mwFreq*10^3);

% Frequency dependence of t and time to frequency mapping
const = trapz(1./v1_range.^2)/t0(end); % const  =  2pi/Qref
t_f = cumtrapz((1/const)*(1./v1_range.^2));
f_adapted = interp1(t_f,f+Params.mwFreq*10^3,t0,'pchip');
f_adapted = f_adapted-Params.mwFreq*10^3;
phi_adapted = 2*pi*cumtrapz(t0,f_adapted); % Phase modulation
phi_adapted = phi_adapted+abs(min(phi_adapted));
A_adapted = A;%interp1(f,A,f_adapted,'pchip');

IQ0 = A.*exp(1i*phi);
IQ0_adapted = A_adapted.*exp(1i*phi_adapted);

err(2) = ~areequal(IQ0,IQ1,1e-12);
err(3) = ~areequal(IQ0_adapted,IQ2,1e-3);


% Sech/tanh and bandwidth-compensated sech/tanh
%--------------------------------------------------------------------------
clear Params Opt
Params.tp = 0.200; % us
Params.TimeStep = 0.00001; % us
Params.Type = 'sech/tanh';
Params.beta = 10; % us
Params.Frequency = [-100 100]; % MHz
Params.Amplitude = 1; % MHz

% Ideal spectrometer magnitude response function
QL = 60; % Q-factor
f0 = 9.2:0.010:9.5; % GHz
dipfreq = 9.35;
v1 = abs(1./(1+1i*QL*(f0/dipfreq-dipfreq./f0)));

Params.mwFreq = 9.34; % GHz
Params.FrequencyResponse = [f0; v1];

[t,IQ1] = pulse(Params);
Opt.Resonator = 'BWcompensation';
[t,IQ2] = pulse(Params,Opt);

% Consistency check in the frequency domain
f = fdaxis(t);
IQ1_ft = abs(fftshift(fft(IQ1)));
IQ1_ft = IQ1_ft/max(IQ1_ft);
IQ2_ft = abs(fftshift(fft(IQ2)));
IQ2_ft = IQ2_ft/max(IQ2_ft);
profile = interp1((f0-Params.mwFreq)*1e3,v1,f);
profile = profile/max(profile);
check = profile.*IQ2_ft;
check = check/max(check);

err(4) = ~areequal(std(IQ1_ft(abs(f)<0.6*Params.Frequency(2))),std(check(abs(f)<0.6*Params.Frequency(2))),1e-2);

% Calculation
t0 = 0:Params.TimeStep:Params.tp;

% Amplitude modulation
A = sech(Params.beta*((t0-Params.tp/2)/Params.tp));
% Frequency modulation
BWinf = (Params.Frequency(2)-Params.Frequency(1))/tanh(Params.beta/2);
f = (BWinf/2)*tanh((Params.beta/Params.tp)*(t0-Params.tp/2));
% Phase modulation
phi = (BWinf/2)*(Params.tp/Params.beta)*log(cosh((Params.beta/Params.tp)*(t0-Params.tp/2)));
phi = 2*pi*phi;

% Bandwidth compensation
v1_range = interp1(f0*10^3,v1,f+Params.mwFreq*10^3);
v1_range = A.*v1_range;

% Frequency dependence of t and time to frequency mapping
const = trapz(f,1./v1_range.^2)/t0(end); % const  =  2pi/Qref
t_f = cumtrapz(f,(1/const)*(1./v1_range.^2));
f_adapted = interp1(t_f,f+Params.mwFreq*10^3,t0,'pchip');
f_adapted = f_adapted-Params.mwFreq*10^3;
phi_adapted = 2*pi*cumtrapz(t0,f_adapted); % Phase modulation
phi_adapted = phi_adapted+abs(min(phi_adapted));
A_adapted = interp1(f,A,f_adapted,'pchip');

IQ0 = A.*exp(1i*phi);
IQ0_adapted = A_adapted.*exp(1i*phi_adapted);

err(5) = ~areequal(IQ0,IQ1,1e-11);
err(6) = ~areequal(IQ0_adapted,IQ2,1e-11);

err = any(err);

data = [];