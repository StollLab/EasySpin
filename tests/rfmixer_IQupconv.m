function [err,data] = test(opt,olddata)

% Check rfmixer() upconversion for IQ data
%--------------------------------------------------------------------------

tp = 0.400; % µs
t = 0:0.0001:tp; % µs
fsweep = 200; % MHz
fcenter = 0; % MHz
k = fsweep/tp;
phi = 2*pi*((fcenter-fsweep/2)*t+(k/2)*t.^2);
signal_bb = exp(1i*phi);

mwFreq = 0.1; % GHz
[tOut,signal_out] = rfmixer(t,signal_bb,mwFreq,'IQshift');

fcenter = mwFreq*1e3;
phi_lo = 2*pi*((fcenter-fsweep/2)*tOut+(k/2)*tOut.^2);
signal_lo = exp(1i*phi_lo);

err = ~areequal(signal_lo,signal_out,1e-12,'abs');

data = [];