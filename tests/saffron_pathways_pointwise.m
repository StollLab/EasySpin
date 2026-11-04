function ok = test() %#ok<FNDEF>

% Indirect dimensions over pulse parameters are computed point by point in
% the subspace engine: echo nutation (flip angle of the first pulse, ideal
% pulses) and a 2D experiment with one delay dimension and one pulse
% dimension, compared with the full-space engine.

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

p90.Flip = pi/2;
p180.Flip = pi;

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0.3 0.2];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p90 0.15 p180 0.15};
Exp.DetPhase = -pi/2;
Exp.Pathways = '+-';

% nutation of the first pulse
E = Exp;
E.nPoints = 9;
E.Dim1 = {'p1.Flip',pi/8};
[x,yP] = saffron(Sys,E,struct('Method','subspace'));
[~,yT] = saffron(Sys,E,struct('Method','fullspace'));
ok(1) = areequal(x,pi/2+pi/8*(0:8),1e-12,'rel');
ok(2) = areequal(yP,yT,1e-3,'rel');

% 2D: phase of the second pulse (point by point) and tau (binned)
E = Exp;
E.nPoints = [4 16];
E.Dim1 = {'p2.Phase',pi/2};
E.Dim2 = {'d1,d2',0.02};
[x,yP] = saffron(Sys,E,struct('Method','subspace','Expand',6));
[~,yT] = saffron(Sys,E,struct('Method','fullspace'));
ok(3) = iscell(x) && isequal(size(yP),[4 16]);
ok(4) = areequal(yP,yT,2e-3,'rel');
