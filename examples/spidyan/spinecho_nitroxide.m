% Running this script will take a few seconds
clear Sys Exp Opt Chirp90 Chirp180

% Spin System
Sys.S = 1/2;
Sys.g = diag([2.00906 2.0062 2.0023]);
Sys.A = diag([11.5 11.5 95]);
Sys.Nucs = '14N';

Exp.Field = 1240; % run the experiment at Q band

% pepper(Sys,Exp) % You can use pepper to see the spectrum and to make sure
% you set your pulses correctly

Symmetry = symm(Sys);
nKnots = 20; % adjust carefully! the more orientations we pick, the better the result, but it might also take much longer

[phi,theta,Weights] = sphgrid(Symmetry,nKnots);

% Pulse echo
% Simulation setup -----------------------------

Chirp90.Type = 'quartersin/linear';
Chirp90.trise = 0.030;
Chirp90.tp = 0.200;
Chirp90.Flip = pi/2;
Chirp90.Frequency = [-0.3 0.3]; % excitation band, GHz

Chirp180.Type = 'quartersin/linear';
Chirp180.trise = 0.030;
Chirp180.tp = 0.100;
Chirp180.Flip = pi;
Chirp180.Frequency = [-0.3 0.3]; % excitation band, GHz

% save time by using five events: the last two events are free evolution
% events, but we only detect during the very last one. This saves time 
% during the simulation (the shorter the detection window, the faster)
Exp.Sequence = {Chirp90 0.25 Chirp180 0.3 0.1}; 

Exp.Flip = [pi/2 pi];
Exp.mwFreq = 34.78; % GHz


% If you want to see only the free evolution after the second pulse, try 
% this instead:
Exp.DetEvents = [0 0 0 0 1];

Opt.DetOperator = {'+1'};
Opt.FreqTranslation = -34.78; % GHz
Opt.SimulationMode = 'LabFrame';


%% A refocused echo with monochromatic pulses

Signal = 0;
for iOrientation = 1 : numel(Weights)
  
  Sys_ = Sys;
  R = erot(phi(iOrientation),theta(iOrientation),0);
  Sys_.g = R*Sys_.g*R';
  Sys_.A = R*Sys_.A*R';
  
  [t, signal_] = spidyan(Sys_,Exp,Opt);
  
  Signal = Signal + signal_*Weights(iOrientation);
  
  progress = [num2str(iOrientation/numel(Weights)*100,'%.2f'),' %'];
  disp(progress)
  
end

%%

figure(2)
clf
plot(t,abs(Signal))
xlabel('t [\mus]')
ylabel('Echo Int. (arb.u.)')