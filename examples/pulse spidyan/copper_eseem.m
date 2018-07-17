% Executing this script will take a couple minutes, the run time can be
% significantly reduced by using a parfor-loop in combination with the 
% parallel computing box

clear

% Spin system
Sys.S = 1/2;
Sys.g = diag([2 2 2.02]); 
Sys.A = diag([20 20 200]); MHz
Sys.Nucs = '63Cu';
Sys.lwpp = 10; 

Exp_.Field = 1240;

% Use this to look at the spectrum (helps with setting the frequency ranges
% of the pulses)
% pepper(Sys,Exp_)

% Get orientations for the loop:
Symmetry = symm(Sys);
nKnots = 40;

[phi,theta,Weights] = sphgrid(Symmetry,nKnots);

% Pulse definitions
Chirp90.Type = 'quartersin/linear';
Chirp90.trise = 0.030;
Chirp90.Frequency = [-0.5 0.5]; % Excitation band, GHz
Chirp90.tp = 0.4;
Chirp90.Flip = pi/2;

Chirp180.Type = 'quartersin/linear';
Chirp180.trise = 0.030;
Chirp180.Frequency = [-0.5 0.5]; % Excitation band, GHz
Chirp180.tp = 0.2;
Chirp180.Flip = pi/2;

Exp.Sequence = {Chirp90 0.25 Chirp180 0.4 0.1}; % Event lengths in us - the lengths of 
                                    % the free evolutions times were chosen
                                    % such, that the echo should be
                                    % centered during the last event
Exp.Field = 1240; % mT
Exp.mwFreq = 35; % Carrier frequency in GHz

% Indirect dimension
Exp.nPoints = 100;
Exp.Dim1 = {'d1,d2', 0.004};

% detection
Exp.DetEvents = [0 0 0 0 1];
Exp.DetOperator = {'+1'};
Exp.DetFreq = Exp.mwFreq; % GHz

%% Orientation loop

Signals = cell(1,numel(Weights));

parfor iOrientation = 1 : numel(Weights)
  
  Sys_ = Sys;
  R = erot(phi(iOrientation),theta(iOrientation),0);
  Sys_.g = R*Sys_.g*R';
  Sys_.A = R*Sys_.A*R';
  
  [Signals{iOrientation}] = spidyan(Sys_,Exp);
  
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

%% Plotting
% Create the time axis
t = 0:(Exp.nPoints-1);
t = t*Exp.Dim{2}+Exp.t(2);

% plot
figure(1)
clf
plot(t,real(Int))
xlabel('\tau (\mus)')
ylabel('Rel. Echo Int. (arb.u.)')