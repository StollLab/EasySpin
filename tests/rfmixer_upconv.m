function ok = test()

% Compare rfmixer() upconversion results with mathematical expression
%--------------------------------------------------------------------------

% With resampling
t = 0:0.1e-3:0.600; % µs
f = 10; % MHz
signalIn = cos(2*pi*f*t);
mwFreq = 0.200000; % GHz

% DSB
[tOut,signal_dsb1] = rfmixer(t,signalIn,mwFreq,'DSB');

signal_dsb2 = 0.5*(cos(2*pi*(f+mwFreq*1e3)*tOut)+cos(2*pi*(f-mwFreq*1e3)*tOut));

% Without resampling

% SSB (upper)
Opt.dt = t(2)-t(1);
[tOut,signal_usb1] = rfmixer(t,signalIn,mwFreq,'USB',Opt);

signal_usb2 = cos(2*pi*(f+mwFreq*1e3)*tOut);

% SSB (lower)
[tOut,signal_lsb1] = rfmixer(t,signalIn,mwFreq,'LSB',Opt);

signal_lsb2 = cos(2*pi*(f-mwFreq*1e3)*tOut);

ok(1) = areequal(signal_dsb1,signal_dsb2,1e-6,'abs');
ok(2) = areequal(signal_usb1,signal_usb2,1e-2,'abs');
ok(3) = areequal(signal_lsb1,signal_lsb2,1e-2,'abs');
