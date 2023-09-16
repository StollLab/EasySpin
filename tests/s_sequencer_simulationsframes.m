function ok = test()

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;

Opt.SimFreq = 31;

[Events1, Vary1] = runprivate('s_sequencer',Exp,Opt);

Pulse.Frequency = 1000* [-0.100 0.100] + 33500;

Exp.Sequence = {Pulse 0.5 Pulse Pulse};

Exp.mwFreq = 0;

[Events2, Vary2] = runprivate('s_sequencer',Exp,Opt);

ok(1) = areequal(Events1{1}.IQ,Events2{1}.IQ,1e-9,'abs');
ok(2) = areequal(Events1{3}.IQ,Events2{3}.IQ,1e-9,'abs');
