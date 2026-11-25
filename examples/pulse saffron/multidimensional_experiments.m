% Setting up indirect dimensions (saffron)
%==========================================================================
% This script shows how to set up a pulse sequence with several pulses and
% then explains on three examples how to run multidimensional experiments

clear

% Spin System
Sys.S = 1/2;

% Experiment
Pulse90.Type = 'rectangular';
Pulse90.tp = 0.05; % µs
Pulse90.Flip = pi/2; % rad

Pulse180.Type = 'rectangular';
Pulse180.tp = 0.05; % µs
Pulse180.Flip = pi; % rad

% pi/2 - tau - pi - tau
Exp.Sequence = {Pulse90 0.5 Pulse180 0.5}; % µs

Exp.mwFreq = 9.5; % GHz
Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.MolFrame = [0 0 0]; % single orientation, molecular frame = lab frame

Exp.DetWindow = 'sequence'; % detect during the whole sequence
Exp.DetOperator = {'z1'};

% Run simulation --------------
[TimeAxis, Signal] = saffron(Sys,Exp);

% Plotting of the initial pulse sequence --------------------
figure(1)
clf
plot(TimeAxis*1000,real(Signal));
xlabel('t (ns)')
ylabel('<S_z>')
axis tight
ylim([-1 1])

%% Two Dimensional Experiment - Increment the flip angle of the second pulse
% The following code changes the flip angle of the second pulse by pi/4
% Exp.nPoints times
Exp.nPoints = 5;
Exp.Dim1 = {'p2.Flip', pi/4}; % rad

[x, Signal] = saffron(Sys,Exp);

% x is a cell array with the axes: the indirect dimension and the time
% axis of the transient; Signal has one row per point of the indirect
% dimension
TimeAxis = x{end};
figure(2)
clf
hold on
for i = 1:size(Signal,1)
  plot(TimeAxis*1000,real(Signal(i,:)));
end
xlabel('t (ns)')
ylabel('<S_z>')
axis tight
ylim([-1 1])

%% Changing the inter pulse delay with a nonlinear increment
% Increment the first delay by 100 and 250 ns. Since this changes the
% length of the detected transient, it can't be done with an indirect
% dimension. Instead, run one simulation per delay.
Exp = rmfield(Exp,{'nPoints','Dim1'});
increments = [0 0.1 0.25]; % µs

figure(3)
clf
hold on
for i = 1:numel(increments)
  Exp_ = Exp;
  Exp_.Sequence{2} = Exp.Sequence{2} + increments(i);
  [TimeAxis, Signal] = saffron(Sys,Exp_);
  plot(TimeAxis*1000,real(Signal));
end
xlabel('t (ns)')
ylabel('<S_z>')
axis tight
ylim([-1 1])

%% Three-dimensional experiment
Exp.nPoints = [2 3]; % 2 steps in 1st and 3 in 2nd dimension
Exp.Dim1 = {'p1.Phase,p2.Phase', pi/4}; % Changes the Phase of both pulses
                                        % by pi/4 each step
Exp.Dim2 = {'p2.Flip', pi/2};  % Change the flip angle of the second pulse

% In order to see the influence of the phase, the detection is changed to
% S^+, detected in the frame rotating at Exp.DetFreq (default Exp.mwFreq)
Exp.DetOperator = {'+1'};
Exp.DetFreq = 9.5; % GHz

[x, Signal] = saffron(Sys,Exp);

% Plotting --------------
TimeAxis = x{end};
figure(4)
clf
hold on
sizeSignal = size(Signal);
LinearSignal = reshape(Signal,[prod(sizeSignal(1:end-1)) sizeSignal(end)]);
for i = 1:size(LinearSignal,1)
  plot(TimeAxis*1000,real(LinearSignal(i,:)));
end
xlabel('t (ns)')
ylabel('<S_x>')
axis tight
