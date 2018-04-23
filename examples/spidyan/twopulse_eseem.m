% A simple ESEEM time trace simulation - running the script will take 
% a while

clear Sys Opt Exp Pulse

% Set up a Gaussian Distribution of (electron) spins
CenterFrequency = 33.5; % center frequency of Gaussian distribution, GHz
GWidth = 0.01; % width of Gaussian distribution, GHz
FreqStart = 33.45;  % starting value for sampling
FreqEnd = 33.55;  % final value for sampling
Sampling = 0.0005;   % stepsize for sampling
ZeemanFreqVec = FreqStart:Sampling:FreqEnd; % vector with resonance frequencies
P = exp(-((CenterFrequency-ZeemanFreqVec)/GWidth).^2); % probabilities
P = P/trapz(P); % normalization

nSpinpackets = length(ZeemanFreqVec);

% Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500;
% Add a nucleus
Sys = nucspinadd(Sys,'1H',[8 45 45]);

%% ESEEM sequence with 
Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.005; % us

Exp.t = [0.05 0.25 0.025 0.125 0.3];
% The first and third event are pulses, the second event is the delay that
% corresponds to the indirect dimension. The fourth event is stepped
% together with the first delay. This way the echo always appears at the
% same position during the last event

Exp.Pulses = {Pulse 0 Pulse 0 0};
Exp.Field = 1240; 
Exp.TimeStep = 0.00001;
Exp.Frequency = [-0.080 0.080];
Exp.Flip = [pi/2 pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = [0 0 0 0 1];

Exp.nPoints = 50;
Exp.Dim = {'d1,d2', 0.004};

Opt.DetOperator = {'+1'};
Opt.FreqTranslation = [-33.5];
Opt.FrameShift = 32;

% Loop over the spinpackets and sum up the traces 
for i = 1 : nSpinpackets
  
  % Replace Sys.ZeemanFreq with ZeemanFreq of current spin packet
  Sys.ZeemanFreq = ZeemanFreqVec(i);
  
  [TimeAxis, Signal] = spidyan(Sys,Exp,Opt);
  
  if i == 1
    TotalSignal = Signal*P(i);
  else
    TotalSignal = TotalSignal + Signal*P(i);
  end
  
  disp([num2str(round(i/nSpinpackets*100,1)) ' %'])
end

%% Data Plotting
% Get maximum of echo at each acquistion point
EchoModulation = zeros(1,Exp.nPoints);
for i = 1 : size(TotalSignal,1)
  EchoModulation(i) = max(abs(TotalSignal(i,:,:)));
end
EchoModulation = EchoModulation - mean(EchoModulation);
EchoModulation = EchoModulation/max(abs(EchoModulation));

% Calculate Time axis of the ESEEM experiment
tau = Exp.t(2)+linspace(0,Exp.Dim{1,2}*(Exp.nPoints-1),Exp.nPoints);

% Plotting of time traces
figure(1)
clf
hold on
plot(tau*1000,EchoModulation)
xlabel('\tau [ns]')
ylabel('Normalized Modulation of Echo Amplitude [a.u.]')
