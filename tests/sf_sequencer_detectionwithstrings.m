function ok = test() %#ok<FNDEF>

% Exp.DetWindow = 'sequence' detects all events of the sequence, without
% appending any, the same as without Exp.DetWindow

Pulse.Flip = pi;
Pulse.tp = 0.03;

Exp.Sequence = {Pulse 0.5 Pulse 0.2};
Exp.Field = 1240;
Exp.mwFreq = 33.5;

Opt = struct;

Exp1 = Exp;
Exp1.DetWindow = 'sequence';
[Events1,~,Opt1] = runprivate('sf_buildevents',Exp1,Opt);

[Events2,~,Opt2] = runprivate('sf_buildevents',Exp,Opt);

det = cellfun(@(E)E.Detection,Events1);
ok(1) = numel(Events1)==numel(Exp.Sequence) && all(det);
ok(2) = ~Opt1.SinglePointDetection;
ok(3) = isequal(Events1,Events2) && isequal(Opt1.SinglePointDetection,Opt2.SinglePointDetection);
