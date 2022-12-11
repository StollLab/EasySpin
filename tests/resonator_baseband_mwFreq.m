function ok = test()

% Check for agreement using resonator() with pulses at baseband or with 
% pulses already at the microwave frequency
%--------------------------------------------------------------------------

% Gaussian pulse at resonant frequency
%--------------------------------------------------------------------------
% Transfer function
ResonatorFrequency = 9.50; % GHz
ResonatorQL = 500;

% Pulse at baseband
%--------------------------------------------------------------------------
Par.tp = 0.100; % us
Par.Flip = pi;
Par.Type = 'gaussian';
Par.tFWHM = 0.050; % us
Par.Frequency = 0;
[tpulse,IQ] = pulse(Par);

mwFreq = 9.50;

% Simulate effect of resonator
[tsim,IQsim] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Resonator compensation
[tcomp,IQcomp] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Pulse at the microwave frequency
%--------------------------------------------------------------------------
Par.Frequency = Par.Frequency+mwFreq*1e3;
[tpulsemw,IQmw] = pulse(Par);

% Simulate effect of resonator
[tsimmw,IQsimmw] = resonator(tpulsemw,real(IQmw),0,ResonatorFrequency,ResonatorQL,'simulate');

% Resonator compensation
[tcompmw,IQcompmw] = resonator(tpulsemw,real(IQmw),0,ResonatorFrequency,ResonatorQL,'compensate');

% Downconvert data at mwFreq
[tsimmw_,IQsimmw_] = rfmixer(tsimmw,IQsimmw,mwFreq,'IQdemod');

% Reposition and phase
splineinterp = @(dt) interp1(tsimmw_+dt,IQsimmw_,tsim,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x(1))*exp(-1i*x(2))-IQsim).^2)/numel(IQsim));
x0 = [0 mean(angle(IQsimmw_))];
x = fminsearch(fitfunc,x0);
IQsimmw_interp = splineinterp(x(1))*exp(-1i*x(2));

ok(1) = areequal(IQsim,IQsimmw_interp,0.05*max(IQsim),'abs');

% Downconvert data at mwFreq
[tcompmw_,IQcompmw_] = rfmixer(tcompmw,IQcompmw,mwFreq,'IQdemod');

% Reposition and phase
splineinterp = @(dt) interp1(tcompmw_+dt,IQcompmw_,tcomp,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x(1))*exp(-1i*x(2))-IQcomp).^2)/numel(IQcomp));
x0 = [-0.0005 mean(angle(IQcompmw_))];
x = fminsearch(fitfunc,x0);
IQcompmw_interp = splineinterp(x(1))*exp(-1i*x(2));

ok(2) = areequal(IQcomp,IQcompmw_interp,0.1,'rel');

% Sech/tanh at resonant frequency
%--------------------------------------------------------------------------
% Transfer function
ResonatorFrequency = 9.50; % GHz
ResonatorQL = 500;

% Pulse at baseband
%--------------------------------------------------------------------------
Par.tp = 0.100; % us
Par.Type = 'sech/tanh';
Par.beta = 10;
Par.Frequency = [-50 50];
[tpulse,IQ] = pulse(Par);

mwFreq = 9.50;

% Simulate effect of resonator
[tsim,IQsim] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Resonator compensation
[tcomp,IQcomp] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Pulse at the microwave frequency
%--------------------------------------------------------------------------
Par.Frequency = Par.Frequency+mwFreq*1e3;
[tpulsemw,IQmw] = pulse(Par);

% Simulate effect of resonator
[tsimmw,IQsimmw] = resonator(tpulsemw,real(IQmw),0,ResonatorFrequency,ResonatorQL,'simulate');

% Resonator compensation
[tcompmw,IQcompmw] = resonator(tpulsemw,real(IQmw),0,ResonatorFrequency,ResonatorQL,'compensate');

% Downconvert data at mwFreq
[tsimmw_,IQsimmw_] = rfmixer(tsimmw,IQsimmw,mwFreq,'IQdemod');

% Reposition and phase
splineinterp = @(dt) interp1(tsimmw_+dt,IQsimmw_,tsim,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x(1))*exp(-1i*x(2))-IQsim).^2)/numel(IQsim));
x0 = [0 0.1];
x = fminsearch(fitfunc,x0);
IQsimmw_interp = splineinterp(x(1))*exp(-1i*x(2));

ok(3) = areequal(IQsim,IQsimmw_interp,0.05*max(IQsim),'abs');

% Downconvert data at mwFreq
[tcompmw_,IQcompmw_] = rfmixer(tcompmw,IQcompmw,mwFreq,'IQdemod');

% Reposition and phase
splineinterp = @(dt) interp1(tcompmw_+dt,IQcompmw_,tcomp,'spline',0);
fitfunc = @(x) sqrt(sum((splineinterp(x(1))*exp(-1i*x(2))-IQcomp).^2)/numel(IQcomp));
x0 = [0 -0.96];
x = fminsearch(fitfunc,x0);
IQcompmw_interp = splineinterp(x(1))*exp(-1i*x(2));

ok(4) = areequal(IQcomp,IQcompmw_interp,0.1*max(IQcomp),'abs');

% Quartersin off-resonance
%--------------------------------------------------------------------------
% Transfer function
ResonatorFrequency = 9.45; % GHz
ResonatorQL = 500;

% Pulse at baseband
%--------------------------------------------------------------------------
Par.tp = 0.100; % us
Par.Type = 'quartersin';
Par.trise = 0.010;
Par.Frequency = 0;
[tpulse,IQ] = pulse(Par);

mwFreq = 9.50;

% Simulate effect of resonator
[tsim,IQsim] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'simulate');

% Resonator compensation
[tcomp,IQcomp] = resonator(tpulse,IQ,mwFreq,ResonatorFrequency,ResonatorQL,'compensate');

% Pulse at the microwave frequency
%--------------------------------------------------------------------------
Par.Frequency = Par.Frequency+mwFreq*1e3;
[tpulsemw,IQmw] = pulse(Par);

% Simulate effect of resonator
[tsimmw,IQsimmw] = resonator(tpulsemw,real(IQmw),0,ResonatorFrequency,ResonatorQL,'simulate');

% Resonator compensation
[tcompmw,IQcompmw] = resonator(tpulsemw,real(IQmw),0,ResonatorFrequency,ResonatorQL,'compensate');

% Downconvert data at mwFreq
[tsimmw_,IQsimmw_] = rfmixer(tsimmw,IQsimmw,mwFreq,'IQdemod');

% Reposition and phase
splineinterp = @(dt) interp1(tsimmw_+dt,IQsimmw_,tsim,'spline');
fitfunc = @(x) sqrt(sum((splineinterp(x(1))*exp(-1i*x(2))-IQsim).^2)/numel(IQsim));
x0 = [0.005 0];
x = fminsearch(fitfunc,x0);
IQsimmw_interp = splineinterp(x(1))*exp(-1i*x(2));

ok(5) = areequal(IQsim,IQsimmw_interp,0.05*max(abs(IQsim)),'abs');

% Downconvert data at mwFreq
[tcompmw_,IQcompmw_] = rfmixer(tcompmw,IQcompmw,mwFreq,'IQdemod');

% Reposition and phase
splineinterp = @(dt) interp1(tcompmw_+dt,IQcompmw_,tcomp,'spline');
fitfunc = @(x) sqrt(sum((splineinterp(x(1))*exp(-1i*x(2))-IQcomp).^2)/numel(IQcomp));
x0 = [0 0];
x = fminsearch(fitfunc,x0);
IQcompmw_interp = splineinterp(x(1))*exp(-1i*x(2));

ok(6) = areequal(IQcomp,IQcompmw_interp,0.1*max(abs(IQcomp)),'abs');
