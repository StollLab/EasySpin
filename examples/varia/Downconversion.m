% Digital downconversion and quadrature detection for a real input signal
% at an intermediate frequency (IF)
%==========================================================================
% In this example, a real signal at an intermediate frequency (IF) is
% downconverted and the quadrature signal is recovered.

clear, clf

% Define signal at an intermediate frequency
%--------------------------------------------------------------------
dt = 0.1e-3; % µs
t = 0:dt:0.500; % µs
x0 = t(end)/2;
w = 0.100; % µs
mwFreq = 0.5; % GHz
echo = gaussian(t,x0,w).*cos(2*pi*mwFreq*1e3*t);

% Digital downconversion and quadrature detection with rfmixer()
%--------------------------------------------------------------------
rfmixer(t,echo,mwFreq,'IQdemod');
