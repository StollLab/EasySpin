% setting up indirect dimensions (spidyan)
%==========================================================================
% This script shows how the set up a pulse sequence with several pulses and
% then explains on four examples how to run multidimensional experiments

clear

% Spin System
Sys.ZeemanFreq = 9.5;

% Experiment
Pulse90.Type = 'rectangular';
Pulse90.tp = 0.05; % µs
Pulse90.Flip = pi; % rad

Pulse180.Type = 'rectangular';
Pulse180.tp = 0.05; % µs
Pulse180.Flip = pi; % rad
  
% pi - tau - pi/2 - tau
Exp.Sequence = {Pulse90 0.5 Pulse180 0.5}; % µs

Exp.mwFreq = 9.5; % GHz

Exp.DetOperator = {'z1'};

% Run simulation --------------
[TimeAxis, Signal] = spidyan(Sys,Exp);

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

[TimeAxis, Signal] = spidyan(Sys,Exp);

% Signal now is an n-dimensional array
figure(2)
clf
hold on
for i = 1:size(Signal,1)
  plot(TimeAxis*1000,real(squeeze(Signal(i,:))));
end
xlabel('t (ns)')
ylabel('<S_z>')
axis tight
ylim([-1 1])
%% Two Dimensional Experiment - Changing the inter pulse delay with a nonlinear increment
% Increment the first delay by 100 and 250 ns
Exp.nPoints = 3; 
Exp.Dim1 = {'d1', [0 0.1 0.25]}; % µs

[TimeAxis, Signal] = spidyan(Sys,Exp);
% Since the length of a detected event is changed, the traces of each
% acquistion point have different lengths and Signal becomes a cell array

% Plotting --------
figure(3)
clf
hold on
for i = 1:numel(Signal)
  % Indexing now corresponds to the elements in the cell array Signal
  plot(TimeAxis{i}*1000,real(Signal{i}));
end
xlabel('t [ns]')
ylabel('<S_z>')
axis tight
ylim([-1 1])

%% Three-dimensional experiment with constant length
Exp.nPoints = [2 3]; % 2 steps in 1st and 3 in 2nd dimension
Exp.Dim1 = {'p1.Phase,p2.Phase', pi/4}; % Changes the Phase of both pulses 
                                        % by pi/4 each step 
Exp.Dim2 = {'p2.Flip', pi/2};  % Change the flip angle of the second pulse

% In order to see the influence of the phase, the detection is changed to
% S^+ and the detected trace down converted with Opt.FreqTranslation
Exp.DetOperator = {'+1'};
Exp.DetFreq = 9.5; % GHz

[TimeAxis, Signal] = spidyan(Sys,Exp);

% Plotting -------------- 
figure(4)
clf
hold on
sizeSignal = size(Signal);
LinearSignal = reshape(Signal,[prod(sizeSignal(1:end-1)) sizeSignal(end)]);
for i = 1:size(LinearSignal,1)
  plot(TimeAxis*1000,real(LinearSignal(i,:)));
end
xlabel('t [ns]')
ylabel('<S_x>')
axis tight

%% 2 Dimensional Experiment with varying length
Exp.nPoints = [2 3]; % 2 steps in 1st and 3 in 2nd dimension
Exp.Dim1 = {'p1.Phase,p2.Phase', pi/4}; % Changes the Phase of both pulses 
                                        % by pi/4 each step 
Exp.Dim2 = {'p2.tp', 0.05};  % Change the length of the second pulse, µs

% Now, Sz and S^+ are detected
Exp.DetOperator = {'z1','+1'};
Exp.DetFreq = [0 9.5]; % GHz

[TimeAxis, Signal] = spidyan(Sys,Exp);

% Plottting -------------- 
figure(5)
clf
hold on
for i = 1:numel(Signal)
  plot(TimeAxis{i}*1000,real(Signal{i}));
end
xlabel('t [ns]')
ylabel('<S_x>')
axis tight
ylim([-1 1])
