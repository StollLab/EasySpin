function ok = test() %#ok<FNDEF>

% Mims ENDOR, single orientation, short rectangular RF pulse with a given
% RF field: the subspace engine (two-level population transfer per NMR
% transition) and the full-space engine (time-stepped RF pulse with the
% nuclear Zeeman operator) agree, including the side lobes of the RF
% excitation profile.

Sys.Nucs = '1H';
Sys.A = [1 1 3];  % MHz
Sys.g = 2.0023;

p90.Flip = pi/2;
rf.Channel = 'rf';
rf.tp = 2;  % µs
rf.B1 = 0.5/(2*0.04258);  % mT, pi pulse for a bare proton

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];
Exp.CrystalSymmetry = 1;
Exp.rfFreq = 13;  % MHz
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
Exp.nPoints = 21;
Exp.Dim1 = {'rfFreq',0.2};
Exp.Pathways = '+00-';
Exp.OnOff = 3;  % ENDOR effect: echo with RF pulse minus echo without
Exp.DetPhase = -pi/2;

[x,y1,info1] = saffron(Sys,Exp,struct('Method','subspace'));
[~,y2,info2] = saffron(Sys,Exp,struct('Method','fullspace'));

ok(1) = strcmp(info1.Method,'subspace') && strcmp(info2.Method,'fullspace');
ok(2) = areequal(real(y1),real(y2),0.03,'rel');
ok(3) = areequal(x,13+0.2*(0:20));
