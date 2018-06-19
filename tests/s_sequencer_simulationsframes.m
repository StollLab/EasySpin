function [err,data] = test(opt,olddata)

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.mwFreq = 33.5;

Opt.SimFrequency = 31;

[Events1, Vary1] = runprivate('s_sequencer',Exp,Opt);

Pulse.Frequency = [-0.100 0.100] + 33.5;

Exp.Sequence = {Pulse 0.5 Pulse Pulse};

Exp = rmfield(Exp,'mwFreq');

[Events2, Vary2] = runprivate('s_sequencer',Exp,Opt);


if any([~areequal(Events1{1}.IQ,Events2{1}.IQ,1e-9) ~areequal(Events1{3}.IQ,Events2{3}.IQ,1e-9)])
  err = 1;
else
  err = 0;
end

data = [];

