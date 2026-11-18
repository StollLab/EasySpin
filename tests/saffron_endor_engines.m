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
Exp.DetPhase = -pi/2;

% echo without RF field (B1 = 0): the same at all RF frequencies, so a
% single point suffices
Eoff = Exp;
Eoff.Sequence{5}.B1 = 0;
Eoff.nPoints = 1;

Opt1 = struct('Method','subspace');
Opt2 = struct('Method','fullspace');
[x,y1,info1] = saffron(Sys,Exp,Opt1);
[~,y2,info2] = saffron(Sys,Exp,Opt2);
[~,y1off] = saffron(Sys,Eoff,Opt1);
[~,y2off] = saffron(Sys,Eoff,Opt2);

% ENDOR effect: reduction of the echo by the RF pulse
endor1 = real(y1off-y1);
endor2 = real(y2off-y2);

ok(1) = strcmp(info1.Method,'subspace') && strcmp(info2.Method,'fullspace');
ok(2) = areequal(endor1,endor2,0.03,'rel') && max(endor1)>0.1*abs(y1off);
ok(3) = areequal(x,13+0.2*(0:20));
ok(4) = areequal(y1off,y2off,1e-4,'rel');

% the sign of the RF field does not change the population transfer
Eneg = Exp;
Eneg.Sequence{5}.B1 = -rf.B1;
[~,y1neg] = saffron(Sys,Eneg,Opt1);
ok(5) = areequal(y1neg,y1,1e-10,'rel');
