% spin echo of a nitroxide using spidyan (spidyan)
%==========================================================================
% this demonstrates how to simulate the spin echo of a nitroxide.
% the same result can be obtained with saffron (see example
% thyme_2p_echo.m)
clear

% Spin System
Sys.S = 1/2;
Sys.g = diag([2.00906 2.0062 2.0023]); % MHz
Sys.A = diag([11.5 11.5 95]); % MHz
Sys.Nucs = '14N';

Exp.Field = 324.9; % mT
% pepper(Sys,Exp) % You can use pepper to see the spectrum and to make sure
% you set your pulses correctly

nKnots = 20; 

Symmetry = symm(Sys);
grid = sphgrid(Symmetry,nKnots);
nOrientations = numel(grid.weights);

% eperiment setup
Chirp90.Type = 'quartersin/linear'; 
Chirp90.trise = 0.030; % smoothed edges, us
Chirp90.tp = 0.200; % pulse length, us
Chirp90.Flip = pi/2; % flip angle in rad
Chirp90.Frequency = [-120 120]; % excitation band, GHz
Chirp90.Phase = pi; 

Chirp180.Type = 'quartersin/linear';
Chirp180.trise = 0.030;% smoothed edges, us
Chirp180.tp = 0.100; % pulse length, us
Chirp180.Flip = pi; % flip angle in rad
Chirp180.Frequency = [-120 120]; % excitation band, GHz
Chirp180.Phase = pi;

% save time by using five events: the last two events are free evolution
% events, but we only detect during the very last one. This saves time 
% during the simulation (the shorter the detection window, the faster)
tau = 0.5; % us
Exp.Sequence = {Chirp90 tau Chirp180 tau+Chirp180.tp}; 
Exp.mwFreq = 9.1; % GHz

Exp.DetWindow = [-0.05 0.05]; % us
Exp.DetOperator = {'+1'};
Exp.DetFreq = 9.1; % GHz


Signal = 0; % initialize Signal
for iOrientation = 1 : nOrientations
  
  Sys_ = Sys; % create temporary copy of Sys
  
  R = erot(grid.phi(iOrientation),grid.theta(iOrientation),0); % rotation matrix
  
  Sys_.g = R'*Sys_.g*R; % rotate Sys.g
  Sys_.A = R'*Sys_.A*R; % rotate Sys.A
  
  [t, signal_] = spidyan(Sys_,Exp);
  
  Signal = Signal + signal_*grid.weights(iOrientation); % sum up signals
  
  progress = [num2str(iOrientation/nOrientations*100,'%.2f'),' %'];
  disp(progress)
  
end

% plotting
figure(1)
clf
plot(t,real(Signal))
xlabel('transient (\mus)')
ylabel('Signal (arb.u.)')
