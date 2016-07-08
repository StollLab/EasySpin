function [err,data] = test(opt,olddata)

% Check rfmixer() upconversion for IQ data
%--------------------------------------------------------------------------

tp = 400; % ns
t = 0:0.1:tp; % ns
fsweep = 0.200; % GHz
fcenter = 0; % GHz
k = fsweep/tp;
phi = 2*pi*((fcenter-fsweep/2)*t+(k/2)*t.^2);
signal_bb = exp(1i*phi);

LOfreq = 0.1; % GHz
[tOut,signal_out_re,signal_out_im] = rfmixer(t,real(signal_bb),imag(signal_bb),LOfreq,'+');
signal_out = signal_out_re + 1i*signal_out_im;

fcenter = LOfreq;
phi_lo = 2*pi*((fcenter-fsweep/2)*tOut+(k/2)*tOut.^2);
signal_lo = exp(1i*phi_lo);

err = ~areequal(signal_lo,signal_out,1e-12);

data = [];