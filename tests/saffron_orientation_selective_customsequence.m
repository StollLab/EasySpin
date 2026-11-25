function ok = test() %#ok<FNDEF>

% Orientation selection with ideal pulses (Opt.IdealPulses, excitation
% window 1/tp): predefined and user-defined two-pulse ESEEM agree

dt = 0.01;
tau = 0;
nPoints = 120;
Field = 350;
mwFreq = 9.5;
tp = 0.005;  % µs, excitation window 200 MHz

Exp.Sequence = '2pESEEM';
Exp.mwFreq = mwFreq;
Exp.dt = dt;
Exp.nPoints = nPoints;
Exp.tp = tp;
Exp.tau = tau;
Exp.Field = Field;

Sys.Nucs = '1H';
Sys.A = [1 6];
Sys.g = [2 2.2];

Opt.GridSize = 91;
Opt.IdealPulses = true;

[~,y1] = saffron(Sys,Exp,Opt);

p90.Flip = pi/2;
p90.tp = tp;
p180.Flip = pi;
p180.tp = tp;

ExpManual.Sequence = {p90 tau p180 tau};
ExpManual.nPoints = nPoints;
ExpManual.Field = Field;
ExpManual.Dim1 = {'d1,d2' dt};
ExpManual.mwFreq = mwFreq;
ExpManual.Pathways = '+-';
ExpManual.DetPhase = -pi/2;

[~, y2] = saffron(Sys,ExpManual,Opt);



ok = areequal(y1,y2,1e-4,'abs');

