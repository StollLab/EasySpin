function [ok,data] = test(opt,olddata)

% Spin system
Sys.S = 1/2;
Sys.g = 1.898468237763169;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

% Experiment
Exp.Sequence = '3pESEEM';
Exp.dt = 0.01;  % µs
Exp.tau = 0.1;  % µs
Exp.T = 0.5;  % µs
Exp.nPoints = 10;
Exp.Field = 350;  % mT

Exp.DetWindow = [-0.05 0.05];  % µs
Exp.mwFreq = 9.3;  % GHz
Exp.DetPhase = 0;

Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 -pi/2 0];
[x1, y1] = saffron(Sys,Exp);

data.x1 = x1;
data.y1 = y1;

Exp.SampleFrame = [0 -pi/2 0; 0 0 0];
[x2, y2] = saffron(Sys,Exp);

data.x2 = x2;
data.y2 = y2;

if ~isempty(olddata)
  ok(1) = areequal(y1,olddata.y1,1e-3,'abs');
  ok(2) = areequal(y2,olddata.y2,1e-3,'abs');
else
  ok = [];
end

if opt.Display
  subplot(2,1,1);
  plot(x1{2},real(y1),x1{2},real(olddata.y1),'.')';
  subplot(2,1,2);
  plot(x2{2},real(y2),x2{2},real(olddata.y2),'.')';
end
