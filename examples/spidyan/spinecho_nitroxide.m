% Running this script will take a few seconds
clear Sys Exp Opt

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

Chirp.Type = 'quartersin/linear';
Chirp.trise = 0.030;


% save time by using five events: the last two events are free evolution
% events, but we only detect during the very last one. This saves time 
% during the simulation (the shorter the detection window, the faster)
Exp.t = [0.200 0.25 0.100 0.3 0.1]; % us 
Exp.Pulses = {Chirp 0 Chirp 0 0}; 

% Since in this simulation we decided to simulate in the rotating frame, we
% need a much shorter timestep! (compare Opt.FrameShift
Exp.TimeStep = 0.00001; % us

Exp.Flip = [pi/2 pi];
Exp.mwFreq = 34.78; % GHz
Exp.Frequency = [-0.3 0.3]; % excitation band, GHz

% If you want to see only the free evolution after the second pulse, try 
% this instead:
Exp.DetEvents = [0 0 0 0 1];

Opt.DetOperator = {'+1'};
Opt.FreqTranslation = -34.78; % GHz


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
plot(t,real(Signal))
xlabel('t [\mus]')
ylabel('Echo Int. (arb.u.)')