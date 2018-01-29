function [err,data] = test(opt,olddata)

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Exp.t = [0.1 0.5 0.1 0.2];
Exp.Pulses = {Pulse 0 Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi pi];
Exp.mwFreq = 33.5;

Opt.FrameShift = 32;

Exp.mwFreq = 0;
Opt = rmfield(Opt,'FrameShift');

[Events1, Vary1] = runprivate('s_sequencer',Exp,Opt);

Exp = rmfield(Exp,'mwFreq');

[Events2, Vary2] = runprivate('s_sequencer',Exp,Opt);


if any([~isequal(Events1,Events2) ~isequal(Vary1,Vary2)])
  err = 1;
else
  err = 0;
end

data = [];

