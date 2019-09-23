function [err,data] = test(opt,olddata)

% the original data sets to compare to where created with the old saffron

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

Field = 324.9;
dt = 0.050;
nPoints = [256 256];
tau = 0.08;

t1 = 0.1;
t2 = 0.1;

Opt.nKnots = 20;

%--- Tests

Exp.Sequence = 'HYSCORE';
Exp.dt = dt;
Exp.nPoints = nPoints;
Exp.tau = tau;
Exp.Field = Field;
Exp.t1 = t1;
Exp.t2 = t2;

[x1, y1] = saffron(Sys,Exp,Opt);

data.x1 = x1;
data.y1 = y1;

% now the same with a manual definition
% The old code, that was used to generate the data:
% ExpManual.Field = Field;
% ExpManual.Flip = [1 1 2 1];
% ExpManual.Inc = [0 1 2 0];
% ExpManual.t = [tau t1 t2 tau];
% ExpManual.dt = dt;
% ExpManual.nPoints = nPoints; 

P90.Flip = pi/2;
P180.Flip = pi;

ExpManual.Sequence = {P90 tau P90 t1 P180 t2 P90 tau};
ExpManual.Field = Field;

ExpManual.nPoints = nPoints;
ExpManual.Dim1 = {'d2' dt};
ExpManual.Dim2 = {'d3' dt};

[x2, y2] = saffron(Sys,ExpManual,Opt);

data.x2 = x2;
data.y2 = y2;

if ~isempty(olddata)
  err = ~areequal(x1{1},olddata.x1,1e-4,'abs') || ...
        ~areequal(y1,olddata.y1,1e-4,'abs') || ...
        ~areequal(x2{1},olddata.x2,1e-4,'abs') || ...
		~areequal(y2,olddata.y2,1e-4,'abs');
else
  err = [];
end

