function [err,data] = test(opt,olddata)

% the original data sets to compare to where created with the old saffron

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

dt = 0.01;
tau = 0;
nPoints = 120;
Field = 350;

Opt.nKnots = 20;

%--- Tests

Exp.Sequence = '2pESEEM';
Exp.dt = 0.01;
Exp.tau = tau;
Exp.nPoints = nPoints;
Exp.Field = Field;

[x1, y1] = saffron2(Sys,Exp,Opt);

data.x1 = x1;
data.y1 = y1;

% now the same with a manual definition
% The old code, that was used to generate the data:
% ExpManual.Field = Field;
% ExpManual.Flip = [1 2];
% ExpManual.Inc = [1 1];
% ExpManual.t = [0 0];
% ExpManual.dt = dt;
% ExpManual.nPoints = nPoints; 

p90.Flip = pi/2;
p180.Flip = pi;

ExpManual.Sequence = {p90 tau p180 tau};
ExpManual.nPoints = nPoints;
ExpManual.Field = Field;
ExpManual.Dim1 = {'d1,d2' dt};

[x2, y2] = saffron2(Sys,ExpManual,Opt);

data.x2 = x2;
data.y2 = y2;

if ~isempty(olddata)
  err = any([~areequal(x1,olddata.x1,1e-4) ~areequal(y1,olddata.y1,1e-4)...
    ~areequal(x2,olddata.x2,1e-4) ~areequal(y2,olddata.y2,1e-4)]);
else
  err = [];
end

