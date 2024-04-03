function [ok,data] = test(opt,olddata)

% Compare pre-defind and user-defined 2-pulse ESEEM

Sys.S = 1/2;
Sys.Nucs = '1H';
Sys.A_ = [5 2];

tau0 = 0;
dt = 0.010;  % µs
nPoints = 120;
Field = 350;  % mT

Opt.GridSize = 20;

ExpPredefined.Sequence = '2pESEEM';
ExpPredefined.dt = dt;  % µs
ExpPredefined.tau = tau0;

ExpPredefined.nPoints = nPoints;
ExpPredefined.Field = Field;

[x1, y1] = saffron(Sys,ExpPredefined,Opt);

p90.Flip = pi/2;
p180.Flip = pi;
ExpManual.Sequence = {p90 tau0 p180 tau0};
ExpManual.Dim1 = {'d1,d2', dt};

ExpManual.nPoints = nPoints;
ExpManual.Field = Field;

[x2, y2] = saffron(Sys,ExpManual,Opt);

data.x1 = x1;
data.y1 = y1;
data.x2 = x2;
data.y2 = y2;

if ~isempty(olddata)
  ok = ...
    areequal(x1,olddata.x1,1e-4,'abs') && ...
    areequal(y1,olddata.y1,1e-4,'abs') && ...
    areequal(x2,olddata.x2,1e-4,'abs') && ...
    areequal(y2,olddata.y2,1e-4,'abs');
else
  ok = [];
end
