% Running this script will take a couple minutes, the processing time can be
% significantly improved by using the parallel computing box
clear Sys Exp Opt

% Spin System
Sys.S = 1/2;
Sys.g = diag([2 2 2.02]);
Sys.A = diag([20 20 200]);
Sys.Nucs = '63Cu';
Sys.lwpp = 10;
Exp.Field = 1240;

% pepper(Sys,Exp)

Symmetry = symm(Sys);
nKnots = 100;

[phi,theta,Weights] = sphgrid(Symmetry,nKnots);

% Pulse echo
% Simulation setup -----------------------------

Chirp.Type = 'quartersin/linear';
Chirp.trise = 0.030;
Exp.Frequency = [-0.5 0.5]; % excitation band, GHz
Exp.t = [0.400 0.25 0.200 0.4 0.1]; % us
Exp.Pulses = {Chirp 0 Chirp 0 0}; 

% Rect.Type = 'rectangular';
% Exp.Frequency = [0]; % GHz
% Exp.t = [0.025 0.25 0.025 0.18 0.2]; % us
% Exp.Pulses = {Rect 0 Rect 0 0}; 

Exp.Field = 1240; % mT
Exp.TimeStep = 0.00005; % us

Exp.Flip = [pi/2 pi];
Exp.mwFreq = 35; % GHz
% If you want to see only the free evolution after the second pulse, try 
% this instead:
Exp.DetEvents = [0 0 0 0 1];

Exp.nPoints = 100;
Exp.Dim = {'d1,d2', 0.004};

Opt.DetOperator = {'+1'};
Opt.FreqTranslation = -35; % GHz
Opt.FrameShift = 30; % Using a frame shift allows us to propagate faster

%% A refocused echo with monochromatic pulses
Signals = cell(1,numel(Weights));

% if you do now have the parallel computing toolbox available use a
% regular for loop:
% for  iOrientation = 1 : numel(Weights)
parfor iOrientation = 1 : numel(Weights)
  
  LoopSys = Sys;
  R = erot(phi(iOrientation),theta(iOrientation),0);
  LoopSys.g = R*LoopSys.g*R';
  LoopSys.A = R*LoopSys.A*R';
  
  [Signals{iOrientation}] = spidyan(LoopSys,Exp,Opt);
  
end

%% Signals need to be summed up
% The par for loop does not allow to cummulate them within the loop

Signal = Signals{1}*Weights(1);

for iOrientation = 2 : numel(Weights)
  Signal = Signals{iOrientation}*Weights(iOrientation);
end

% find the maximum of each transient
Int = zeros(1,Exp.nPoints);
for iTrace = 1 : Exp.nPoints
  Int(iTrace) = max(abs(Signal(iTrace,1,:)));
end

Int = Int/max(Int);

% Create the time axis
t = 0:(Exp.nPoints-1);
t = t*Exp.Dim{2}+Exp.t(2);


% plot
figure(1)
clf
plot(t,real(Int))
xlabel('\tau [\mus]')
ylabel('Rel. Echo Int. a.u.')