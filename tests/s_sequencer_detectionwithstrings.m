function ok = test()

Pulse.Flip = pi;
Pulse.tp = 0.03;

Exp.Sequence = {Pulse 0.5 Pulse 0.2};
Exp.Field = 1240; 

Exp.mwFreq = 33.5;

Opt.SimulationMode = 'step wise';

% First Test - All Events -------------------
Exp1 = Exp;
Exp1.DetSequence = 1;
[Events1] = runprivate('s_sequencer',Exp1,Opt);

Exp2 = Exp;
Exp2.DetSequence = 'all';
[Events2] = runprivate('s_sequencer',Exp2,Opt);

% Second Test - Single Point after a pulse -------------------

Exp3 = Exp;
Exp3.DetSequence = [0 0 0 1];

[Events3] = runprivate('s_sequencer',Exp3,Opt);

Exp4 = Exp;
Exp4.DetSequence = 'last';
[Events4] = runprivate('s_sequencer',Exp4,Opt);

ok = all([isequal(Events1,Events2) isequal(Events3,Events4)]);
