function [ok,data] = test(opt,olddata)

% the original data sets to compare to where created with the old saffron

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

dt = 0.01;
tau = 0.1;
nPoints = 120;
Field = 350;

Opt.GridSize = 20;

%--- Tests

Exp.Sequence = '3pESEEM';
Exp.dt = dt;
Exp.tau = tau;
Exp.nPoints = nPoints;
Exp.Field = Field;
Exp.MolFrame = [0 0 0];

Exp.SampleFrame = [0 -pi/2 0];

[x1, y1] = saffron(Sys,Exp,Opt);

data.x1 = x1;
data.y1 = y1;

Exp.SampleFrame = [0 -pi/2 0; 0 0 0];

[x2, y2] = saffron(Sys,Exp,Opt);

data.x2 = x2;
data.y2 = y2;

if ~isempty(olddata)
  ok = areequal(y1,olddata.y1,1e-4,'abs') && areequal(y2,olddata.y2,1e-4,'abs');
else
  ok = [];
end

