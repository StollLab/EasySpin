clear Exp Sys Opt Pulse
% See the effect of a (hypothetical) phase cycle. 
% First run the script as it is, and it will give you an inversion pulse
% followed by a free evolution period.
% Then add the phasecycle by uncommenting the line starting with 
% Exp.PhaseCycle.
% Since the detection phase in the provided phase cycle are 1 and -1, the
% sum of the two phasecycles adds up to zero.
% Phase Cycles can be assigned to every pulse in your sequence.
% Connected phasecycles (phase cycling two pulses together) are not
% possible with Exp.PhaseCycle, but can be written by adding an additional
% dimension (see documentation) where the phase of pulses is changed.
% Merging along this dimension has to be done manually.

% Spin System
Sys.S = 1/2;
Sys.ZeemanFreq = 33.500; % GHz

% Pulse Definitions
Pulse.Type = 'rectangular';

% Define the phase cycle with phase and detection phase 
PhaseCycle = [0   1; 
              pi -1];

% Experiment/Sequence
Exp.t = [0.025 0.4]; % us
Exp.Frequency = 0; % GHz
Exp.Pulses = {Pulse 0};
Exp.Flip = pi;
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Add the phasecycle to the experiment structure
Exp.PhaseCycle = {PhaseCycle};

% Options
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32; % GHz

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting 
figure(1)
clf
plot(TimeAxis,Signal)
ylim([-1 1])
xlabel('t [ns]')
ylabel('<S_z>')