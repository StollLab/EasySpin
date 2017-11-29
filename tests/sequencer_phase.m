function [err,data] = test(opt,olddata)

Pulse.Type = 'quartersin/linear';
Pulse.trise = 0.015; % us

Pulse1  = Pulse;
Pulse1.Phase = pi;

Pulse2  = Pulse;
Pulse2.Phase = pi/2;

Exp.t = [0.1 0.5 0.1 0.2];
Exp.Pulses = {Pulse1 0 Pulse Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = [-0.100 0.100];
Exp.Flip = [pi pi pi];
Exp.mwFreq = 33.5;

Opt.FrameShift = 32;

% First Test ---------------------------------------

[Events1] =  runprivate2('s_sequencer',Exp,Opt);

Exp.Phase = pi;
Exp.Pulses = {Pulse 0 Pulse Pulse};

[Events2] = runprivate2('s_sequencer',Exp,Opt);

% Second Test ---------------------------------------

Exp = rmfield(Exp,'Phase');
Exp.Pulses = {Pulse1 0 Pulse1 Pulse};

[Events3] = runprivate2('s_sequencer',Exp,Opt);

Exp.Phase = [pi/2 pi];
Exp.Pulses = {Pulse2 0 Pulse Pulse};

[Events4] = runprivate2('s_sequencer',Exp,Opt);

% ---------------------------------------------------

if any([~isequal(Events1,Events2) ~isequal(Events3,Events4)])
  err = 1;
else
  err = 0;
end

data = [];

