function ok = test()

Pulse.Flip = pi;
Pulse.tp = 0.03;

Exp.Sequence = {Pulse 0.5 Pulse 0.2};
Exp.Field = 1240; 

Exp.mwFreq = 33.5;

% First Test - Single Point at end of free evolution -------------------
Exp1 = Exp;
Exp1.DetWindow = 0;
Opt = struct;
[Events1] = runprivate('s_sequencer',Exp1,Opt);

Exp2 = Exp;
Exp2.DetWindow = [0 0];
[Events2] = runprivate('s_sequencer',Exp2,Opt);

Exp3 = Exp;
Exp3.Sequence = {Pulse 0.5 Pulse 0.2 0};
Exp3.DetSequence = [0 0 0 0 1];
[Events3] = runprivate('s_sequencer',Exp3,Opt);

% Second Test - Single Point after a pulse -------------------

Exp4 = Exp;
Exp4.Sequence = {Pulse 0.5 Pulse};
Exp4.DetWindow = 0;

[Events4] = runprivate('s_sequencer',Exp4,Opt);

Exp5 = Exp4;
Exp5.DetWindow = [0 0];
[Events5] = runprivate('s_sequencer',Exp5,Opt);

Exp6 = Exp;
Exp6.Sequence = {Pulse 0.5 Pulse 0};
Exp6.DetSequence = [0 0 0 1];

[Events6] = runprivate('s_sequencer',Exp6,Opt);

% Second Test - Single Point somewhere -------------------

% With Pulse as the last event
Exp7 = Exp;
Exp7.Sequence = {Pulse 0.5 Pulse};
Exp7.DetWindow = 0.5;

[Events7] = runprivate('s_sequencer',Exp7,Opt);

Exp8 = Exp;
Exp8.Sequence = {Pulse 0.5 Pulse 0.5 0};
Exp8.DetSequence = [0 0 0 0 1];

[Events8] = runprivate('s_sequencer',Exp8,Opt);

% Free evolution
Exp9 = Exp;
Exp9.DetWindow = -0.1;

[Events9] = runprivate('s_sequencer',Exp9,Opt);

Exp10 = Exp;
Exp10.Sequence = {Pulse 0.5 Pulse 0.1 0};
Exp10.DetSequence = [0 0 0 0 1];

[Events10] = runprivate('s_sequencer',Exp10,Opt);

ok(1) = isequal(Events1,Events2);
ok(2) = isequal(Events1,Events3);
ok(3) = isequal(Events4,Events5);
ok(4) = isequal(Events4,Events6);
ok(5) = isequal(Events7,Events8);
ok(6) = isequal(Events9,Events10);
