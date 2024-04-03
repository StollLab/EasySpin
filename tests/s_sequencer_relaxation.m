function ok = test()

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % µs
Pulse.tp = 0.1;
Pulse.Frequency = 1000* [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.3 Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % µs
Exp.mwFreq = 33.5;
Exp.DetSequence = [1 0 0 0]; 

Opt = struct;

% First Test ---------------------------------------

Opt.Relaxation = 1;

[Events1] = runprivate('s_sequencer',Exp,Opt);

Opt.Relaxation = [1 1 1 1];

[Events2] = runprivate('s_sequencer',Exp,Opt);


ok = isequal(Events1,Events2);
