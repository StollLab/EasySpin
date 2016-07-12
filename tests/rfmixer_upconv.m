function [err,data] = test(opt,olddata)

% Compare rfmixer() upconversion results with mathematical expression
%--------------------------------------------------------------------------

% With resampling
t = 0:1:600; % ns
f = 0.010; % GHz
signalIn = cos(2*pi*f*t);
LOfreq = 0.500; % GHz

% DSB
[tOut,signal_dsb1] = rfmixer(t,signalIn,LOfreq,'DSB');

signal_dsb2 = 0.5*(cos(2*pi*(f+LOfreq)*tOut)+cos(2*pi*(f-LOfreq)*tOut));

% Without resampling
t = 0:0.1:600; % ns
f = 0.010; % GHz
signalIn = cos(2*pi*f*t);
LOfreq = 0.500; % GHz

% SSB (upper)
opt.dt = t(2)-t(1);
[tOut,signal_usb1] = rfmixer(t,signalIn,LOfreq,'SSBup',opt);

signal_usb2 = cos(2*pi*(f+LOfreq)*tOut);

% SSB (lower)
[tOut,signal_lsb1] = rfmixer(t,signalIn,LOfreq,'SSBdown',opt);

signal_lsb2 = cos(2*pi*(f-LOfreq)*tOut);

err(1) = ~areequal(signal_dsb1,signal_dsb2,1e-6);
err(2) = ~areequal(signal_usb1,signal_usb2,1e-2);
err(3) = ~areequal(signal_lsb1,signal_lsb2,1e-2);

err = any(err);

data = [];