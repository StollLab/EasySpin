function ok = test()

Pulse.Flip = pi;
Pulse.tp = 0.03;

Exp.Sequence = {Pulse 0.5 Pulse 0.2};
Exp.Field = 1240; 

Exp.mwFreq = 33.5;

Opt = struct;

% First Test ---------------------------------------
Exp1 = Exp;
Exp1.DetWindow = [0 0.5];

[Events1] = runprivate('s_sequencer',Exp1,Opt);

Exp2 = Exp;
Exp2.Sequence = {Pulse 0.5 Pulse 0.2 0.5};
Exp2.DetSequence = [0 0 0 0 1];

[Events2] = runprivate('s_sequencer',Exp2,Opt);

% Second Test ---------------------------------------
Exp3 = Exp;
Exp3.DetWindow = [-0.20 0.5];

[Events3] = runprivate('s_sequencer',Exp3,Opt);

Exp4 = Exp;
Exp4.Sequence = {Pulse 0.5 Pulse 0 0.7};
Exp4.DetSequence = [0 0 0 0 1];

[Events4] = runprivate('s_sequencer',Exp4,Opt);

% Third Test ---------------------------------------
Exp5 = Exp;
Exp5.Sequence = {Pulse 0.5 Pulse};
Exp5.DetWindow = [0.2 0.5];

[Events5] = runprivate('s_sequencer',Exp5,Opt);

Exp6 = Exp;
Exp6.Sequence = {Pulse 0.5 Pulse 0.2 0.3};
Exp6.DetSequence = [0 0 0 0 1];

[Events6] = runprivate('s_sequencer',Exp6,Opt);

% Fourth Test ---------------------------------------
Exp7 = Exp;
Exp7.Sequence = {Pulse 0.5 Pulse};
Exp7.DetWindow = [0 0.5];

[Events7] = runprivate('s_sequencer',Exp7,Opt);

Exp8 = Exp;
Exp8.Sequence = {Pulse 0.5 Pulse 0.5};
Exp8.DetSequence = [0 0 0 1];

[Events8] = runprivate('s_sequencer',Exp8,Opt);

ok = all([isequal(Events1,Events2) isequal(Events3,Events4) ...
     isequal(Events5,Events6) isequal(Events7,Events8)]);
