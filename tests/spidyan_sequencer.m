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

Opt.ComplexExcitation = [1 0];
Opt.Relaxation = [1 0];
Opt.StateTrajectories = 1;

[Events1] = sequencer(Exp,Opt);

Opt.ComplexExcitation = [1 0 0];
Opt.Relaxation = [1 0 0];
Exp.DetEvents = [1 0 0]; 
Opt.StateTrajectories = [1 1 1 1];

[Events2] = sequencer(Exp,Opt);

% Second Test ---------------------------------------

Exp = rmfield(Exp,'DetEvents');
Opt.StateTrajectories = [1 1];

[Events3] = sequencer(Exp,Opt);

Exp.DetEvents = 0; 
Opt.StateTrajectories = [1 1 0];

[Events4] = sequencer(Exp,Opt);

% Third Test ---------------------------------------

Exp.nPoints = [3 3];
Exp.Dim1 = {'d1', 0.2};
Exp.Dim2 = {'p1.t', 0.2};

[~, Vary1] = sequencer(Exp,Opt);

Exp.Dim1 = {'d1', [0.2 0.4]};
Exp.Dim2 = {'p1.t', [0.2 0.4]};

[~, Vary2] = sequencer(Exp,Opt);

% Fourth Test ---------------------------------------
Exp.t(end) = 0.8;
Exp.Pulses = {Pulse 0 Pulse};
Exp.nPoints = 3;
Exp.Dim1 = {'p2.Position', 0.2};

[~, Vary3] = sequencer(Exp,Opt);

Exp.Dim1 = {'p2.Position', [0.2 0.4]};

[~, Vary4] = sequencer(Exp,Opt);

% Fifth Test ---------------------------------------
Exp.mwFreq = 0;
Opt = rmfield(Opt,'FrameShift');

[Events5, Vary5] = sequencer(Exp,Opt);

Exp = rmfield(Exp,'mwFreq');

[Events6, Vary6] = sequencer(Exp,Opt);


if any([~isequal(Events1,Events2) ~isequal(Events3,Events4) ...
    ~isequal(Vary1,Vary2) ~isequal(Vary3,Vary4) ...
    ~isequal(Events5,Events6) ~isequal(Vary5,Vary6)])
  err = 1;
else
  err = 0;
end

data = [];

