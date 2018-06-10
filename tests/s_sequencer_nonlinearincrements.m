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
Exp.DetEvents = [1 0 0 0]; 

Opt = [];

% First Test ---------------------------------------

Exp.nPoints = [3 3];
Exp.Dim1 = {'d1', 0.2};
Exp.Dim2 = {'p1.t', 0.2};

[~, Vary1] = runprivate('s_sequencer',Exp,Opt);

Exp.Dim1 = {'d1', [0 0.2 0.4]};
Exp.Dim2 = {'p1.t', [0 0.2 0.4]};

[~, Vary2] = runprivate('s_sequencer',Exp,Opt);

% Second Test ---------------------------------------

Exp.t(end) = 0.8;
Exp.Pulses = {Pulse 0 Pulse 0};
Exp.nPoints = 3;
Exp.Dim1 = {'p2.Position', 0.2};

[~, Vary3] = runprivate('s_sequencer',Exp,Opt);

Exp.Dim1 = {'p2.Position', [0 0.2 0.4]};

[~, Vary4] = runprivate('s_sequencer',Exp,Opt);


if any([~isequal(Vary1,Vary2) ~isequal(Vary3,Vary4)])
  err = 1;
else
  err = 0;
end

data = [];

