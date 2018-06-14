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

Opt.FrameShift = 31;
Opt.SimulationMode = 'FrameShift';

Exp.mwFreq = 0;
Opt = rmfield(Opt,'FrameShift');
Opt = rmfield(Opt,'SimulationMode');

[Events1, Vary1] = runprivate('s_sequencer',Exp,Opt);

Exp = rmfield(Exp,'mwFreq');

[Events2, Vary2] = runprivate('s_sequencer',Exp,Opt);


if any([~isequal(Events1,Events2) ~isequal(Vary1,Vary2)])
  err = 1;
else
  err = 0;
end

data = [];

