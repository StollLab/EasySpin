% Upconversion of an IQ signal to an intermediate frequency (IF)
%==========================================================================
% In this example, a sech/tanh pulse defined with in-phase and quadrature
% data is upconverted to an intermediate frequency.

clear, clf

% Define IQ input signal (sech/tanh pulse)
%--------------------------------------------------------------------
Par.tp = 0.200; % µs
Par.Type = 'sech/tanh'; % sech/tanh pulse
Par.beta = 10.6; % truncation parameter
Par.Frequency = [-100 100]; % MHz 
[t,y] = pulse(Par);

% Upconversion to an intermediate frequency
%--------------------------------------------------------------------
mwFreq = 0.100; % GHz
rfmixer(t,y,mwFreq,'IQshift');