% Running this script will take a few seconds

clear Exp Sys Opt Pulse

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
Pulse.Type = 'rectangular';

Exp.t = [0.025 0.25 0.05 0.5]; % us
Exp.Pulses = {Pulse 0 Pulse 0}; 
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [0]; % GHz
Exp.Flip = [pi/2 pi];
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

figure(1)
clf
plot(t*1000,abs(Signal));
xlabel('t [ns]')
axis tight
ylim([0 1])

%% Echo with linear chirps

% Experiment setup -------------------------
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.005; % us

Exp.t = [0.05 0.25 0.025 0.5];
Exp.Frequency = [-0.080 0.080];
Exp.Pulses = {Pulse 0 Pulse 0};

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