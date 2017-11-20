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
Exp.DetEvents = [1 0]; 

Opt.FrameShift = 32;

% First Test ---------------------------------------

Opt.StateTrajectories = 1;

[Events1] = sequencer(Exp,Opt);

Opt.StateTrajectories = [1 1 1 1];

[Events2] = sequencer(Exp,Opt);

% Second Test ---------------------------------------

Opt.StateTrajectories = [1 0];

[Events3] = sequencer(Exp,Opt);

Opt.StateTrajectories = [1 0 0 0];

[Events4] = sequencer(Exp,Opt);


if any([~isequal(Events1,Events2) ~isequal(Events3,Events4)])
  err = 1;
else
  err = 0;
end

data = [];

