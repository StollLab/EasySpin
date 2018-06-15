function [err,data] = test(opt,olddata)

Pulse.Flip = pi;
Pulse.tp = 0.03;

Exp.Sequence = {Pulse 0.5 Pulse 0.2};
Exp.Field = 1240; 

Exp.mwFreq = 33.5;



% First Test ---------------------------------------
Exp1 = Exp;
Exp1.DetWindow = [0 0.5];

[Events1] = runprivate('s_sequencer',Exp1);

Exp2 = Exp;
Exp2.Sequence = {Pulse 0.5 Pulse 0.2 0.5};
Exp2.DetSequence = [0 0 0 0 1];

[Events2] = runprivate('s_sequencer',Exp2);

% Second Test ---------------------------------------
Exp3 = Exp;
Exp3.DetWindow = [-0.20 0.5];

[Events3] = runprivate('s_sequencer',Exp3);

Exp4 = Exp;
Exp4.Sequence = {Pulse 0.5 Pulse 0 0.7};
Exp4.DetSequence = [0 0 0 0 1];

[Events4] = runprivate('s_sequencer',Exp4);

% Third Test ---------------------------------------
Exp5 = Exp;
Exp5.Sequence = {Pulse 0.5 Pulse};
Exp5.DetWindow = [0.2 0.5];

[Events5] = runprivate('s_sequencer',Exp5);

Exp6 = Exp;
Exp6.Sequence = {Pulse 0.5 Pulse 0.2 0.3};
Exp6.DetSequence = [0 0 0 0 1];

[Events6] = runprivate('s_sequencer',Exp6);

% Fourth Test ---------------------------------------
Exp7 = Exp;
Exp7.Sequence = {Pulse 0.5 Pulse};
Exp7.DetWindow = [0 0.5];

[Events7] = runprivate('s_sequencer',Exp7);

Exp8 = Exp;
Exp8.Sequence = {Pulse 0.5 Pulse 0.5};
Exp8.DetSequence = [0 0 0 1];

[Events8] = runprivate('s_sequencer',Exp8);

if any([~isequal(Events1,Events2) ~isequal(Events3,Events4) ...
     ~isequal(Events5,Events6) ~isequal(Events7,Events8)])
  err = 1;
else
  err = 0;
end

data = [];

