function [err,data] = test(opt,olddata)

Pulse.Type = 'rectangular';

Exp.t = [0.1 0.5 0.1];
Exp.Pulses = {Pulse 0 Pulse};
Exp.Field = 1240; 
Exp.TimeStep = 0.0001; % us
Exp.Frequency = 0;
Exp.Flip = [pi pi];
Exp.mwFreq = 33.5;
Exp.DetEvents = 1; 

Opt.FrameShift = 32;

% This tests different ways of defining a dimension where the frequency is
% swept - They should all give the same results

% First Syntax -----------------------------
Exp.nPoints = 3;
Exp.Dim1 = {'p1.Frequency' 0.05};

[~, Vary1] = runprivate('s_sequencer',Exp,Opt);

% Second Syntax -----------------------------
Exp.Dim1 = {'p1.Frequency' [0.05 0.05]};

[~, Vary2] = runprivate('s_sequencer',Exp,Opt);

% Third Syntax -----------------------------
Exp.Dim1 = {'p1.Frequency' [0 0; 0.05 0.05; 0.1 0.1]};

[~, Vary3] = runprivate('s_sequencer',Exp,Opt);

% Fourth Syntax -----------------------------
Exp.Frequency = [0 0];
Exp.Dim1 = {'p1.Frequency(1),p1.Frequency(2)' 0.05};

[~, Vary4] = runprivate('s_sequencer',Exp,Opt);

% Fourth Syntax -----------------------------
Exp.Frequency = [0 0];
Exp.Dim1 = {'p1.Frequency(1)' 0.05;
           'p1.Frequency(2)' 0.05};

[~, Vary5] = runprivate('s_sequencer',Exp,Opt);

if any([~isequal(Vary1,Vary2) ~isequal(Vary1,Vary3) ~isequal(Vary1,Vary4) ~isequal(Vary1,Vary5)])
  err = 1;
else
  err = 0;
end

data = [];