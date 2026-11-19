function ok = test() %#ok<FNDEF>

% Detected transients whose length (number of points) changes along an
% indirect dimension are not supported. The propagation stops at the first
% point whose length differs, with an error that names the workaround.
% Dimensions that keep the detected length constant work.

Sys.S = 1/2;

Pulse.Type = 'rectangular';
Pulse.tp = 0.02;  % µs
Pulse.Flip = pi/2;

Exp.Field = mhz2mt(9500);  % mT, resonance field for g = gfree at 9.5 GHz
Exp.mwFreq = 9.5;  % GHz
Exp.MolFrame = [0 0 0];  % single orientation
Exp.DetOperator = {'z1'};
Exp.nPoints = 3;

% Whole-sequence detection with an incremented delay
E = Exp;
E.Sequence = {Pulse 0.1};
E.DetWindow = 'sequence';
E.Dim1 = {'d1' 0.05};
ok(1) = throwsWith(@() saffron(Sys,E),'one simulation per point');

% Whole-sequence detection with an incremented pulse length
E.Dim1 = {'p1.tp' 0.01};
ok(2) = throwsWith(@() saffron(Sys,E),'one simulation per point');

% Whole-sequence detection with an incremented pulse phase: constant length
E.Dim1 = {'p1.Phase' pi/2};
[~,y] = saffron(Sys,E);
ok(3) = isnumeric(y) && size(y,1)==3;

% A detection window has a fixed length
E = Exp;
E.Sequence = {Pulse 0.1};
E.DetWindow = [0 0.05];  % µs
E.Dim1 = {'d1' 0.05};
[~,y] = saffron(Sys,E);
ok(4) = isnumeric(y) && size(y,1)==3;

end

function ok = throwsWith(fcn,str)
try
  fcn();
  ok = false;
catch err
  ok = contains(err.message,str);
end
end
