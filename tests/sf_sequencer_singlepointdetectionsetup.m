function ok = test() %#ok<FNDEF>

% A single number in Exp.DetWindow, or two equal numbers, gives single-point
% detection: a detected zero-length event at the given time after the end
% of the sequence (after a delay or after a pulse)

tol = 1e-12;

Pulse.Flip = pi;
Pulse.tp = 0.03;

Exp.Sequence = {Pulse 0.5 Pulse 0.2};
Exp.Field = 1240;
Exp.mwFreq = 33.5;

Opt = struct;

% Single point at the end of the last delay
Exp1 = Exp;
Exp1.DetWindow = 0;
[Events1,len,det,single] = events(Exp1,Opt);
Exp2 = Exp;
Exp2.DetWindow = [0 0];
Events2 = events(Exp2,Opt);
ok(1) = isequal(Events1,Events2);
ok(2) = single && areequal(len,[-1 0.5 -1 0.2 0],tol,'abs') && isequal(det,[0 0 0 0 1]);

% Single point right after a pulse
Exp4 = Exp;
Exp4.Sequence = {Pulse 0.5 Pulse};
Exp4.DetWindow = 0;
[Events4,len,det,single] = events(Exp4,Opt);
Exp5 = Exp4;
Exp5.DetWindow = [0 0];
Events5 = events(Exp5,Opt);
ok(3) = isequal(Events4,Events5);
ok(4) = single && areequal(len,[-1 0.5 -1 0],tol,'abs') && isequal(det,[0 0 0 1]);

% Single point some time after a pulse
Exp7 = Exp4;
Exp7.DetWindow = 0.5;
[~,len,det,single] = events(Exp7,Opt);
ok(5) = single && areequal(len,[-1 0.5 -1 0.5 0],tol,'abs') && isequal(det,[0 0 0 0 1]);

% Single point during the last delay
Exp9 = Exp;
Exp9.DetWindow = -0.1;
[~,len,det,single] = events(Exp9,Opt);
ok(6) = single && areequal(len,[-1 0.5 -1 0.1 0],tol,'abs') && isequal(det,[0 0 0 0 1]);

end

% Events, lengths of the free evolution events (-1 for pulses), detection
% flags, and whether the detection is a single point
function [Events,len,det,single] = events(Exp,Opt)
[Events,~,Opt] = runprivate('sf_buildevents',Exp,Opt);
len = -ones(1,numel(Events));
det = zeros(1,numel(Events));
for k = 1:numel(Events)
  if strcmp(Events{k}.type,'free evolution')
    len(k) = Events{k}.t(end);
  end
  det(k) = Events{k}.Detection;
end
single = Opt.SinglePointDetection;
end
