function [err,data] = test(opt,olddata)

% This test makes sure, that 'Frequency' is correctly incremented.
% Since the frequencies are currently defined in GHz, and pulse() requires
% MHz, they need to be converted by s_sequencer. There was a bug, where
% Exp.Dim1 = {'pX.Frequency' [0 0.005]} would not be converted to MHz. This
% test is here to make sure, that bug is not reintroduced. Bonus: it should
% also catch if we change the interface from GHz to MHz and forget to adapt
% s_sequencer correctly

TimeStep = 0.0001;
Flip = pi;
InitFrequency  = 1; % in GHz
df = -0.05; % frequency increment
tp = 0.01;

Pulse.Type = 'rectangular';
Pulse.tp = tp
Pulse.Flip = Flip;
Pulse.Frequency = InitFrequency;

% For manually calling pulse --------
Manual.Type = 'rectangular';  
Manual.TimeStep = TimeStep;
Manual.tp = tp;
Manual.Flip = Flip;
% pulse(Manual)
% -------------

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = TimeStep; % us

Opt.SimulationMode = 'LabFrame';

% First Test - Change one pulse along one dimension -----------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Frequency' df};

[~, Vary1] = runprivate('s_sequencer',Exp,Opt);

Seq1IQs = Vary1.Pulses{1}.IQs;

IQs = cell(1,Exp.nPoints);
for i = 1:3
  Manual.Frequency = (InitFrequency + (i-1)*df)*1000;
  [~,IQs{i}] = pulse(Manual);
end

% Second Test - Change both pulses in the same dimension ------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Frequency,p2.Frequency' df};

[~, Vary2] = runprivate('s_sequencer',Exp,Opt);

Seq2IQs1 = Vary2.Pulses{1}.IQs;
Seq2IQs2 = Vary2.Pulses{2}.IQs;

% same input, but different syntax
Exp.Dim1 = {'p1.Frequency' df;
            'p2.Frequency' df};
          
[~, Vary2_] = runprivate('s_sequencer',Exp,Opt);

% Third Test - Two-dimensional --------------------------------------------
Exp.nPoints = [3 3];
Exp.Dim1 = {'p1.Frequency' 2*df};
Exp.Dim2 = {'p1.Frequency' df};

[~, Vary3] = runprivate('s_sequencer',Exp,Opt);

Seq3IQs = Vary3.Pulses{1}.IQs;

IQs3 = cell(Exp.nPoints);
for i = 1:3
  for j = 1:3
  Manual.Frequency = (InitFrequency + (i-1)*2*df + (j-1)*df)*1000;
  [~,IQs3{i,j}] = pulse(Manual);
  end
end

% -------------------------------------------------------------------------

if any([~isequal(IQs,Seq1IQs) ~isequal(IQs,Seq2IQs1) ~isequal(IQs,Seq2IQs2)...
    ~isequal(Vary2,Vary2_) ~isequal(IQs3{3},Seq3IQs{3})])
  err = 1;
else
  err = 0;
end

data = [];

plot = false; 

if plot
  for i = 1:9
    figure(i+10)
    clf
    hold on
    plot(real(IQs3{i}))
    plot(real(Seq3IQs{i}))
  end
end