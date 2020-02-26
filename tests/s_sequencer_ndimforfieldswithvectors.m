function ok = test()

orig_state = warning;
warning('off','all');

Pulse.Type = 'sech/uniformQ';  

Pulse.beta = 10;
Pulse.n = [6 6];
Pulse.Frequency = 1000* [0.95 1.05];
Pulse.Flip = pi;
Pulse.tp = 0.1;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240;
Exp.mwFreq = 0;
Opt.IntTimeStep = 0.0001; % us

Exp.DetSequence = 1;

% First Syntax -----------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.n(2)' -2};

Opt.SimFreq = 0;
Opt.SimulationMode = 'step wise';

[~, Vary1] = runprivate('s_sequencer',Exp,Opt);

Seq1IQs = Vary1.Pulses{1}.IQs;

% Second Syntax -----------------------------
Exp.Dim1 = {'p1.n(2)' [0 -2 -4]};

[~, Vary2] = runprivate('s_sequencer',Exp,Opt);

Seq2IQs = Vary2.Pulses{1}.IQs;

Pulse.n = [6 6];
Pulse.TimeStep = 0.0001;
Pulse.tp = 0.1;
Pulse.Frequency = [0.95 1.05]*1000;
Pulse.Flip = pi;
% pulse(Pulse)

IQs = cell(1,Exp.nPoints);
for i = 1:3
  Pulse.n(2) = 6 - (i-1)*2;
  [~,IQs{i}] = pulse(Pulse);
end

ok(1) = isequal(IQs,Seq1IQs);
ok(2) = isequal(IQs,Seq2IQs);

warning(orig_state);