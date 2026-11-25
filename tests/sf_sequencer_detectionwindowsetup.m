function ok = test() %#ok<FNDEF>

% Exp.DetWindow appends a detected free evolution to the sequence: after a
% delay, the window start shortens or lengthens the last delay; after a
% pulse, a positive window start adds an undetected delay. Windows that
% start in a pulse are rejected (by saffron).

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

% Windows and single points that start in a pulse: rejected by saffron,
% for both engines
Sys.g = 2;
Exp5 = Exp3;
Exp5.DetWindow = [-0.1 0.5];
Exp6 = Exp;
Exp6.DetWindow = [-0.3 0.5];
Exp7 = Exp;
Exp7.DetWindow = -0.3;
ok(5) = throws(@() saffron(Sys,Exp5),'during a pulse') && ...
  throws(@() saffron(Sys,Exp6),'into a pulse') && ...
  throws(@() saffron(Sys,Exp7,struct('Method','manifold')),'into a pulse');

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

function ok = throws(fcn,msg)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,msg);
end
end
