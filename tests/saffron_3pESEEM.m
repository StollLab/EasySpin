function [err,data] = test(opt,olddata)

% the original data sets to compare to where created with the old saffron

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

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

% now the same with a manual definition
% The old code, that was used to generate the data:
% ExpManual.Field = 350;
% ExpManual.Flip = [1 1 1];
% ExpManual.Inc = [0 1 0];
% ExpManual.t = [tau 0 tau];
% ExpManual.dt = dt;
% ExpManual.nPoints = nPoints; 

p90.Flip = pi/2;

ExpManual.Sequence = {p90 tau p90 0 p90 tau};
ExpManual.Field = Field;

ExpManual.nPoints = nPoints;
ExpManual.Dim1 = {'d2' dt};

[x2, y2] = saffron2(Sys,ExpManual,Opt);

data.x2 = x2;
data.y2 = y2;

if ~isempty(olddata)
  err = any([~areequal(x1,olddata.x1,1e-4) ~areequal(y1,olddata.y1,1e-4)...
    ~areequal(x2,olddata.x2,1e-4) ~areequal(y2,olddata.y2,1e-4)]);
else
  err = [];
end

