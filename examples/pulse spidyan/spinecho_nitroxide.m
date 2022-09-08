% Spin echo of a nitroxide using spidyan (spidyan)
%==========================================================================
% This demonstrates how to simulate the spin echo of a nitroxide with
% shapd pulses.

clear

% Define spin system
Sys.S = 1/2;
Sys.g = diag([2.00906 2.0062 2.0023]);
Sys.A = diag([11.5 11.5 95]);  % MHz
Sys.Nucs = '14N';

Exp.Field = 324.9;  % mT
% Use pepper to simulate the frequency-swept spectrum and to make sure
% pulse excitation bands are set appropriately.
pepper(Sys,Exp)

%%
% Set up pulses
Chirp90.Type = 'quartersin/linear'; 
Chirp90.trise = 0.030;  % rise time for smoothed edges, µs
Chirp90.tp = 0.200;  % pulse length, µs
Chirp90.Flip = pi/2;  % flip angle, rad
Chirp90.Frequency = [-120 120];  % excitation band, MHz
Chirp90.Phase = pi; 

Chirp180.Type = 'quartersin/linear';
Chirp180.trise = 0.030;  % rise time for smoothed edges, µs
Chirp180.tp = 0.100;  % pulse length, µs
Chirp180.Flip = pi;  % flip angle, rad
Chirp180.Frequency = [-120 120];  % excitation band, MHz
Chirp180.Phase = pi;

% Define pulse sequence
tau = 0.5;  % delay between pulses, µs
Exp.Sequence = {Chirp90, tau, Chirp180, tau+Chirp180.tp};

Exp.mwFreq = 9.1;  % center excitation frequency, GHz
Exp.DetFreq = 9.1;  % detection frequency, GHz

% Detect transient over time window around last time point in sequence
Exp.DetWindow = [-0.05 0.05]; % start and end time, µs

% Set up orientational grid for powder averaging
GridSize = 20;
Symmetry = symm(Sys);
grid = sphgrid(Symmetry,GridSize);
nOrientations = numel(grid.weights);

Signal = 0;
for iOrientation = 1:nOrientations
  
  
  % Rotate all tensors of spin system
  Sys_ = Sys; % create temporary copy of Sys
  R = erot(grid.phi(iOrientation),grid.theta(iOrientation),0); % rotation matrix
  Sys_.g = R'*Sys_.g*R; % rotate Sys.g
  Sys_.A = R'*Sys_.A*R; % rotate Sys.A
  
  % Simulate signal and accumulate
  [t, signal_] = spidyan(Sys_,Exp);
  Signal = Signal + signal_*grid.weights(iOrientation); % accumulate signals
  
  progress = [num2str(iOrientation/nOrientations*100,'%.2f'),' %'];
  disp(progress)
  
end

% Plotting
figure(1)
clf
plot(t,real(Signal),t,imag(Signal));
xlabel('time relative to sequence start time (µs)');
ylabel('signal (arb.u.)');
legend('in phase (I)','out-of-phase (Q)')
