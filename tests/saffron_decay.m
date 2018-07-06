function [err,data] = test(opt,olddata)

% the original data sets to compare to where created with the old saffron

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];
Sys.T1 = 2;
Sys.T2 = 1;

dt = 0.01;
tau = 0.1;
nPoints = 120;
Field = 350;

Opt.nKnots = 20;

%--- Tests

Exp.Sequence = '3pESEEM';
Exp.dt = dt;
Exp.tau = tau;
Exp.nPoints = nPoints;
Exp.Field = Field;

[x1, y1] = saffron2(Sys,Exp,Opt);

data.x1 = x1;
data.y1 = y1;

if ~isempty(olddata)
  err = any([~areequal(x1,olddata.x1,1e-4) ~areequal(y1,olddata.y1,1e-4)]);
else
  err = [];
end

