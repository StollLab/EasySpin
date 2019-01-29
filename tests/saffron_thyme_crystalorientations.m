function [err,data] = test(opt,olddata)

Sys.S = 1/2;
Sys.g = 1.898468237763169;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

dt = 0.01;
tau = 0.1;
nPoints = 10;
Field = 350;

%--- Tests

Exp.Sequence = '3pESEEM';
Exp.dt = dt;
Exp.tau = tau;
Exp.T = 0.5;
Exp.nPoints = nPoints;
Exp.Field = Field;

Exp.DetWindow = [-0.05 0.05];
Exp.mwFreq = 9.3;
Exp.DetPhase = 0;

Exp.CrystalOrientation = [0 pi/2 0];

[x1, y1] = saffron(Sys,Exp);

data.x1 = x1;
data.y1 = y1;

Exp.CrystalOrientation = [0 pi/2 0; 0 0 0];

[x2, y2] = saffron(Sys,Exp);

data.x2 = x2;
data.y2 = y2;

if ~isempty(olddata)
  err = any([~areequal(y1,olddata.y1,1e-4) ~areequal(y2,olddata.y2,1e-4)]);
else
  err = [];
end

