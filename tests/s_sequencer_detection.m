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

% First Test ---------------------------------------

Opt.DetEvents = true;

[Events1] = runprivate('s_sequencer',Exp,Opt);

Opt.DetEvents = [1 1 1 1];

[Events2] = runprivate('s_sequencer',Exp,Opt);

% Second Test ---------------------------------------

Opt.DetEvents = false;

[Events3] = runprivate('s_sequencer',Exp,Opt);

Opt.DetEvents = [0 0 0 0];

[Events4] = runprivate('s_sequencer',Exp,Opt);

Opt.DetEvents = [];

[Events5] = runprivate('s_sequencer',Exp,Opt);

Opt = rmfield(Opt,'DetEvents');

[Events6] = runprivate('s_sequencer',Exp,Opt);

% Third Test ---------------------------------------

Opt.DetEvents = [1 1];

[Events7] = runprivate('s_sequencer',Exp,Opt);

Opt.DetEvents = [1 1 0 0];

[Events8] = runprivate('s_sequencer',Exp,Opt);


if any([~isequal(Events1,Events2) ~isequal(Events3,Events4) ...
    ~isequal(Events3,Events4) ~isequal(Events3,Events5) ...
    ~isequal(Events3,Events6) ~isequal(Events7,Events8)])
  err = 1;
else
  err = 0;
end

data = [];

