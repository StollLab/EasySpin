function [err,data] = test(opt,olddata)

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us
Pulse.tp = 0.1;
Pulse.Frequency = [-0.100 0.100];
Pulse.Flip = pi;

Exp.Sequence = {Pulse 0.5 Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us;
Exp.mwFreq = 33.5;
Exp.DetEvents = [1 0 0 0]; 

% First Test ---------------------------------------

Opt.StateTrajectories = 1;

[Events1] = runprivate('s_sequencer',Exp,Opt);


Opt.StateTrajectories = [1 1 1 1];

[Events2] = runprivate('s_sequencer',Exp,Opt);


if ~isequal(Events1,Events2)
  err = 1;
else
  err = 0;
end

data = [];

