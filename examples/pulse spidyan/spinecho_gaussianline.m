% Running this script will take a few seconds

clear Exp Sys Opt Pulse90 Pulse 180

% Sets up a Gaussian Distribution of spin packets
CenterFrequency = 33.5; % center frequency of Gaussian distribution, GHz
GWidth = 0.01;     % width of Gaussian distribution, GHz
FreqStart = 33.45;  % starting value for sampling
FreqEnd = 33.55;  % final value for sampling
Sampling = 0.0005;   % stepsize for sampling
ZeemanFreqVec = FreqStart:Sampling:FreqEnd; % vector with resonance frequencies
P = exp(-((CenterFrequency-ZeemanFreqVec)/GWidth).^2); % probabilities
P = P/trapz(P); % normalization
nSpinpackets = length(ZeemanFreqVec);

% Spin System
Sys.S = 1/2;

%% A refocused echo with monochromatic pulses

% Simulation setup -----------------------------
Pulse90.Type = 'rectangular';
Pulse90.tp = 0.025;
Pulse90.Flip = pi/2;

Pulse180.Type = 'rectangular';
Pulse180.tp = 0.025;
Pulse180.Flip = pi/2;

Exp.Sequence = {Pulse90 0.25 Pulse90 0.5}; 
Exp.Field = 1240; % mT
Exp.mwFreq = 33.5; % GHz

% This detects the entire experiment:
Exp.DetEvents = 1;
% If you want to see only the free evolution after the second pulse, try 
% this instead:
Exp.DetEvents = [0 0 0 1];

Opt.DetOperator = {'+1'};
Opt.FrameShift = 32; % GHz
Opt.SimulationMode = 'FrameShift';
Opt.FreqTranslation = -33.5; % GHz


Signal = 0;
% Loop over the spinpackets and sum up the traces ------------
for i = 1 : nSpinpackets
  
  Sys.ZeemanFreq = ZeemanFreqVec(i); % Set Zeeman frequency
  
  [t, signal] = spidyan(Sys,Exp,Opt);
  
  Signal = Signal + signal*P(i);
  
end

figure(1)
clf
plot(t*1000,abs(Signal));
xlabel('t [ns]')
axis tight
ylim([0 1])

%% Echo with linear chirps

% Experiment setup -------------------------
Chirp90.Type = 'quartersin/linear';
Chirp90.trise = 0.005; % us
Chirp90.tp = 0.05; % us
Chirp90.Frequency = [-0.080 0.080];

Chirp180.Type = 'quartersin/linear';
Chirp180.trise = 0.005; % us
Chirp180.tp = 0.025; % us
Chirp180.Frequency = [-0.080 0.080];

Exp.Sequence = {Chirp90 0.25 Chirp180 0.5};

% Loop over the spinpackets and sum up the traces ------------
for i = 1 : nSpinpackets
  
  Sys.ZeemanFreq = ZeemanFreqVec(i); % Set Zeeman frequency
  
  [t, signal] = spidyan(Sys,Exp,Opt);
  if i == 1
    Signal = signal*P(i);
  else
    Signal = Signal + signal*P(i);
  end
  
end

figure(2)
clf
plot(t*1000,abs(Signal));
xlabel('t [ns]')
axis tight
ylim([0 1])