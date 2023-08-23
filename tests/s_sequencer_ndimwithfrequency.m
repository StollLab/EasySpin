function ok = test()

Pulse.Type = 'rectangular';
Pulse.tp = 0.1;
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.DetSequence = 1; 

Opt.SimFreq = 32;
Opt.SimulationMode = 'step wise';

% This tests different ways of defining a dimension where the frequency is
% swept - They should all give the same results

% First Syntax -----------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Frequency' 0.05};

[~, Vary1] = runprivate('s_sequencer',Exp,Opt);

% Second Syntax -----------------------------
Exp.Dim1 = {'p1.Frequency' [0.05 0.05]};

[~, Vary2] = runprivate('s_sequencer',Exp,Opt);

% Third Syntax -----------------------------
Exp.Dim1 = {'p1.Frequency' [0 0; 0.05 0.05; 0.1 0.1]};

[~, Vary3] = runprivate('s_sequencer',Exp,Opt);

% Fourth Syntax -----------------------------
Pulse.Frequency = 1000* [0 0];
Exp.Sequence = {Pulse 0.5 Pulse};

Exp.Dim1 = {'p1.Frequency(1),p1.Frequency(2)' 0.05};

[~, Vary4] = runprivate('s_sequencer',Exp,Opt);

% Fifth Syntax -----------------------------
Exp.Frequency = 1000* [0 0];
Exp.Dim1 = {'p1.Frequency(1)' 0.05;
           'p1.Frequency(2)' 0.05};

[~, Vary5] = runprivate('s_sequencer',Exp,Opt);

ok(1) = isequal(Vary1,Vary2);
ok(2) = isequal(Vary1,Vary3);
ok(3) = isequal(Vary1,Vary4);
ok(4) = isequal(Vary1,Vary5);
