function [err,data] = test(opt,olddata)
% Check amplitude for simulation of transmitter nonlinearity in 
% transmitter()
%--------------------------------------------------------------------------

% Effect of transmitter nonlinearity on pulse shape and amplitude
%--------------------------------------------------------------------------
nu1_max = 15;

% Standard sech amplitude modulation
Params.Type = 'sech/none';
Params.TimeStep = 0.00025; % us
Params.tp = 0.200; % us
Params.beta = 10;
Params.Amplitude = nu1_max;

[t,IQ] = pulse(Params);

% Simulation of the effect of amplitude compression
InputAmplitude = 0.75;
Ain =  0:0.001:1;
Aout = nu1_max*(Ain - 0.30*Ain.^3); % MHz

IQ_compressed = transmitter(InputAmplitude*IQ/max(IQ),Ain,Aout,'simulate');

% Compare output amplitude with amplitude expected based on the defined
% transmitter gain curve
[~,ind] = min(abs(Ain-InputAmplitude));
err = ~areequal(max(abs(IQ_compressed)),Aout(ind),1e-11);

data = [];