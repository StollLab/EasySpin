function ok = test() %#ok<FNDEF>

% Exp.DetPhase with several detection operators: one phase per operator
% gives the same signals as single-operator runs with those phases, a
% single value applies to all operators, a wrong number of values is an
% error, and an empty DetOperator means the default.

Sys.S = 1/2;
Sys.g = gfree*mhz2mt(9505)/mhz2mt(9500);  % 5 MHz above the carrier

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;

Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.mwFreq = 9.5;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.Sequence = {Pulse 0.05};
Exp.DetWindow = 0;

phases = [pi/3 -pi/5];
ops = {'+1','x1'};

% Both operators at once, each with its own phase
E = Exp;
E.DetOperator = ops;
E.DetPhase = phases;
y = saffron(Sys,E);

% Each operator alone
yRef = zeros(1,2);
for k = 1:2
  E = Exp;
  E.DetOperator = ops(k);
  E.DetPhase = phases(k);
  yRef(k) = saffron(Sys,E);
end
ok(1) = areequal(y(:),yRef(:),1e-10,'abs');

% A single value applies to all operators
E = Exp;
E.DetOperator = ops;
E.DetPhase = phases(1);
y1 = saffron(Sys,E);
E.DetPhase = phases([1 1]);
y2 = saffron(Sys,E);
ok(2) = areequal(y1,y2,1e-12,'abs');

% Wrong number of values
E.DetPhase = [0 0 0];
ok(3) = throwsWith(@() saffron(Sys,E),'Exp.DetPhase');

% An empty DetOperator means the default (S+)
E = Exp;
E.DetOperator = {'+1'};
yPlus = saffron(Sys,E);
E.DetOperator = [];
yEmpty = saffron(Sys,E);
ok(4) = areequal(yEmpty,yPlus,1e-10,'abs');

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
