function ok = test() %#ok<FNDEF>

% Exp.DetWindow appends a detected free evolution to the sequence: after a
% delay, the window start shortens or lengthens the last delay; after a
% pulse, a positive window start adds an undetected delay. Windows that
% start in a pulse are rejected.

tol = 1e-12;

Pulse.Flip = pi;
Pulse.tp = 0.03;

Exp.Sequence = {Pulse 0.5 Pulse 0.2};
Exp.Field = 1240;
Exp.mwFreq = 33.5;

Opt = struct;

% Sequence ending with a delay, window starting at its end
Exp1 = Exp;
Exp1.DetWindow = [0 0.5];
[len,det] = events(Exp1,Opt);
ok(1) = areequal(len,[-1 0.5 -1 0.2 0.5],tol,'abs') && isequal(det,[0 0 0 0 1]);

% Window starting before the end of the last delay
Exp2 = Exp;
Exp2.DetWindow = [-0.2 0.5];
[len,det] = events(Exp2,Opt);
ok(2) = areequal(len,[-1 0.5 -1 0 0.7],tol,'abs') && isequal(det,[0 0 0 0 1]);

% Sequence ending with a pulse, window starting after it
Exp3 = Exp;
Exp3.Sequence = {Pulse 0.5 Pulse};
Exp3.DetWindow = [0.2 0.5];
[len,det] = events(Exp3,Opt);
ok(3) = areequal(len,[-1 0.5 -1 0.2 0.3],tol,'abs') && isequal(det,[0 0 0 0 1]);

% Sequence ending with a pulse, window starting at its end
Exp4 = Exp3;
Exp4.DetWindow = [0 0.5];
[len,det] = events(Exp4,Opt);
ok(4) = areequal(len,[-1 0.5 -1 0.5],tol,'abs') && isequal(det,[0 0 0 1]);

% Windows that start in a pulse
Exp5 = Exp3;
Exp5.DetWindow = [-0.1 0.5];
Exp6 = Exp;
Exp6.DetWindow = [-0.3 0.5];
ok(5) = throws(@() runprivate('sf_buildevents',Exp5,Opt)) && ...
  throws(@() runprivate('sf_buildevents',Exp6,Opt));

end

% Lengths of the free evolution events (-1 for pulses) and detection flags
function [len,det] = events(Exp,Opt)
Events = runprivate('sf_buildevents',Exp,Opt);
len = -ones(1,numel(Events));
det = zeros(1,numel(Events));
for k = 1:numel(Events)
  if strcmp(Events{k}.type,'free evolution')
    len(k) = Events{k}.t(end);
  end
  det(k) = Events{k}.Detection;
end
end

function ok = throws(fcn)
try
  fcn();
  ok = false;
catch
  ok = true;
end
end
