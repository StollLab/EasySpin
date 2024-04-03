function ok = test()

dt = 0.01;
tau = 0;
nPoints = 120;
Field = 350;
mwFreq = 9.5;
ExciteWidth = 200;

Exp.Sequence = '2pESEEM';
Exp.mwFreq = mwFreq;
Exp.dt = dt;
Exp.nPoints = nPoints;
Exp.ExciteWidth = ExciteWidth;
Exp.tau = tau;
Exp.Field = Field;

Sys.Nucs = '1H';
Sys.A = [1 6];
Sys.g = [2 2.2];

Opt.GridSize = 91;

[~,y1] = saffron(Sys,Exp,Opt);

p90.Flip = pi/2;
p180.Flip = pi;

ExpManual.Sequence = {p90 tau p180 tau};
ExpManual.nPoints = nPoints;
ExpManual.Field = Field;
ExpManual.Dim1 = {'d1,d2' dt};
ExpManual.mwFreq = mwFreq;
ExpManual.ExciteWidth = ExciteWidth;

[~, y2] = saffron(Sys,ExpManual,Opt);



ok = areequal(y1,y2,1e-4,'abs');

