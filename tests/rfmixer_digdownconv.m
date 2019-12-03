function [err,data] = test(opt,olddata)

% Check rfmixer() digital downconversion
%--------------------------------------------------------------------------

dt = 0.0002; % µs
t = 0:dt:0.500; % µs
x0 = t(end)/2;
w = 0.100; % µs
mwFreq = 0.5; % GHz
amplitude = gaussian(t,x0,w);
modulation = cos(2*pi*mwFreq*1e3*t);

echo = amplitude.*modulation;

[tOut,signalOut] = rfmixer(t,echo,mwFreq,'IQdemod');

err = ~areequal(amplitude,real(signalOut),1e-7,'abs');

data = [];