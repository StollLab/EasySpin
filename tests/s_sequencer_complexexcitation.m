function ok = test()

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.Flip = pi;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.tp = 0.1;

Pulse2 = Pulse;
Pulse2.tp = 0.2;

Exp.Sequence = {Pulse 0.5 Pulse Pulse2};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % µs

Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0 0 0]; 

Opt.SimulationMode = 'step wise';

% First Test ---------------------------------------

[Events1] = runprivate('s_sequencer',Exp,Opt);

Exp.mwPolarization = 'circular';

[Events2] = runprivate('s_sequencer',Exp,Opt);


ok = ~isequal(Events1,Events2);
