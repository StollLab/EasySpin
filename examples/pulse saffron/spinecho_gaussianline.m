% Spin echo of a Gaussian line (saffron)
%==========================================================================
% this demonstrates how to simulate the spin echo of a Gaussian line, as a
% sum over spin packets with different resonance frequencies, each one a
% single orientation. This is useful when you are only interested in the
% dynamics during your pulse sequence and don't need to simulate a powder
% spectrum, and it can speed up simulations.
% This script calculates two echos:
%  - first with rectangular pulses
%  - second with frequency-swept pulses.

clear

%% set up a Gaussian line of spin packets
CenterFrequency = 9.5; % center frequency of Gaussian distribution, GHz
GWidth = 0.02;     % width of Gaussian distribution, GHz
FreqStart = 9.45;  % starting value for sampling
FreqEnd = 9.55;  % final value for sampling
Sampling = 0.0010;   % stepsize for sampling
ZeemanFreqVec = FreqStart:Sampling:FreqEnd; % vector with resonance frequencies
P = exp(-((CenterFrequency-ZeemanFreqVec)/GWidth).^2); % probabilities
P = P/trapz(P); % normalization
nSpinpackets = length(ZeemanFreqVec);

%% refocused echo with monochromatic pulses

% Experiment definition
Pulse90.Type = 'rectangular';
Pulse90.tp = 0.025; % pulse length, µs
Pulse90.Flip = pi/2; % flip angle, rad

Pulse180.Type = 'rectangular';
Pulse180.tp = 0.025; % pulse length, µs
Pulse180.Flip = pi;  % flip angle, rad

Exp.Sequence = {Pulse90 0.1 Pulse180};
Exp.mwFreq = 9.5; % GHz
Exp.Field = mhz2mt(CenterFrequency*1e3);  % mT, resonance field of the line center for g = gfree
Exp.MolFrame = [0 0 0]; % single orientation, molecular frame = lab frame

% Detect the free evolution during 200 ns after the second pulse. To
% observe the entire sequence, append the delay to Exp.Sequence and use
% Exp.DetWindow = 'sequence'.
Exp.DetWindow = [0 0.2]; % µs, relative to the end of the sequence

Exp.DetOperator = {'+1'}; % detect electron coherence
Exp.DetFreq = 9.5; % GHz
% Exp.DetPhase = pi; % rad, for proper phasing of the signal

Signal = 0; % initialize the signal

% Loop over the spinpackets and sum up the traces ------------
for i = 1 : nSpinpackets

  Sys_.g = gfree*mhz2mt(ZeemanFreqVec(i)*1e3)/Exp.Field;  % g of the spin packet

  [t, signal] = saffron(Sys_,Exp);

  Signal = Signal + signal*P(i); % cumulate signals

end

% plotting
figure
subplot(2,1,1)
plot(t,real(Signal),t,imag(Signal));
xlabel('time after the last pulse (\mus)')
axis tight


%% echo with linear chirps

% Experiment definition
Chirp90.Type = 'quartersin/linear';
Chirp90.trise = 0.005; % rise time in µs
Chirp90.tp = 0.05; % pulse length, µs
Chirp90.Frequency = [-80 80]; % frequency band
Chirp90.Flip = pi/2; % flip angle in rad

Chirp180.Type = 'quartersin/linear';
Chirp180.trise = 0.005; % rise time in µs
Chirp180.tp = 0.025; % pulse length, µs
Chirp180.Frequency = [-80 80]; % frequency band
Chirp180.Flip = pi; % flip angle in rad

Exp.Sequence = {Chirp90 0.1 Chirp180}; % free evolution events in µs
% Exp.DetPhase = 0; % rad, for proper phasing of the signal

Signal = 0; % initialize the signal

% Loop over the spinpackets and sum up the traces
for i = 1 : nSpinpackets

  Sys_.g = gfree*mhz2mt(ZeemanFreqVec(i)*1e3)/Exp.Field;  % g of the spin packet

  [t, signal] = saffron(Sys_,Exp);
  Signal = Signal + signal*P(i);

end

subplot(2,1,2)
plot(t,real(Signal),t,imag(Signal));
xlabel('time after the last pulse (\mus)')
axis tight
