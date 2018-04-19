clear Exp Sys Opt Pulse
% This script shows how the set up a pulse sequence with several pulses and
% then explains on four examples how to run multidimensional experiments

% Spin System -----------------------
Sys.S = [1/2];
Sys.ZeemanFreq = [33.500];

Pulse.Type = 'rectangular';

% This Exp structure creates a pi - tau - pi/2 - tau pulse sequence ------
Exp.t = [0.05 0.5 0.05 0.5]; % us
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; % mT
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0; % GHz
Exp.Flip = [pi pi/2];
Exp.mwFreq = 33.5; % GHz
Exp.DetEvents = 1;

% Options -------------------------
Opt.DetOperator = {'z1'};
Opt.FrameShift = 32;

% Run simulation --------------
[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plotting of the initial pulse sequence --------------------
figure(1)
clf
plot(TimeAxis*1000,real(Signal));
xlabel('t [ns]')
ylabel('<S_z>')
axis tight
ylim([-1 1])

%% 1 Dimensional Experiment - Increment the flip angle of the second pulse
% The following code changes the flip angle of the second pulse by pi/4
% Exp.nPoints times
Exp.nPoints = 5;
Exp.Dim1 = {'p2.Flip', pi/4};

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);
% Signal now is a three dimensional array, with the first dimension being
% the index of the acquisition (nPoints), second the index of the detection 
% operator and third the detected points


% Plotting -----------------------
figure(2)
clf
hold on
for i = 1:size(Signal,1)
  plot(TimeAxis*1000,real(squeeze(Signal(i,:,:))));
end
xlabel('t [ns]')
ylabel('<S_z>')
axis tight
ylim([-1 1])


%% 1 Dimensional Experiment - Changing the inter pulse delay with a nonlinear increment
% Increment the first delay by 100 and 250 ns
Exp.nPoints = 3; 
Exp.Dim1 = {'d1', [0.1 0.25]};  % The first data point is always the one 
                                % defined in the Exp structure, hence only 
                                % [0.1 0.25] is written, even if 
                                % Exp.nPoints = 3

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);
% Since the length of a detected event is changed, the traces of each
% acquistion point have different lengths and Signal has to be a cell array

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


%% 2 Dimensional Experiment with constant length
Exp.nPoints = [2 3]; % 2 steps in 1st and 3 in 2nd dimension
Exp.Dim1 = {'p1.Phase,p2.Phase', pi/4}; % Changes the Phase of both pulses 
                                        % by pi/4 each step 
Exp.Dim2 = {'p2.Flip', pi/2};  % Change the flip angle of the second pulse

% In order to see the influence of the phase, the detection is changed to
% S^+ and the detected trace down converted with Opt.FreqTranslation
Opt.DetOperator = {'+1'};
Opt.FreqTranslation = [-33.5]; % GHz

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

% Plottting -------------- 
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
Exp.Dim2 = {'p2.tp', 0.05};  % Change the length of the second pulse

% Now, Sz and S^+ are detected
Opt.DetOperator = {'z1','+1'};
Opt.FreqTranslation = [0 -33.5]; % GHz

[TimeAxis, Signal] = spidyan(Sys,Exp,Opt);

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
