function [err,data] = test(opt,olddata)

% Check rfmixer() digital downconversion
%--------------------------------------------------------------------------

dt = 0.2; % ns
t = 0:dt:500; % ns
x0 = t(end)/2;
w = 100; % ns
IFfreq = 0.5; % GHz
amplitude = gaussian(t,x0,w);
modulation = cos(2*pi*IFfreq*t);

echo = amplitude.*modulation;

[tOut,signalOut] = rfmixer(t,echo,IFfreq,'IQdemod');

err = ~areequal(amplitude,real(signalOut));

data = [];