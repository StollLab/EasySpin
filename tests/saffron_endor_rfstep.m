function ok = test() %#ok<FNDEF>

% Full engine, RF frequency sweep: the time step of the RF waveform
% follows the RF frequency of each point. A sweep that starts at 1 MHz and
% jumps to 12 MHz agrees on the common points with a linear sweep starting
% at 12 MHz. The RF pulse has its full length, also when 40 steps per RF
% period do not fit into it a whole number of times. The RF pulse at 1 MHz
% is far from the NMR transitions, so that echo is the echo without RF
% pulse, and the ENDOR effects are compared relative to it.

Sys.Nucs = '1H';
Sys.A = [1 1 3];  % MHz

p90.Flip = pi/2;
rf.Channel = 'rf';
rf.tp = 1.33;  % µs, 638.4 RF periods/40 at 12 MHz
rf.B2 = 0.5/(rf.tp*0.04258);  % mT, pi pulse for a bare proton

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0 0];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p90 0.2 p90 0 rf 0.1 p90 0.2};
Exp.Pathways = '+00-';
Exp.DetPhase = -pi/2;

Opt.Method = 'full';

% 1 MHz, then 12 ... 16 MHz
ExpA = Exp;
ExpA.rfFreq = 1;  % MHz
ExpA.Dim1 = {'rfFreq',[0 11:0.5:15]};
ExpA.nPoints = 10;
[xA,yA] = saffron(Sys,ExpA,Opt);

% 12 ... 16 MHz
ExpB = Exp;
ExpB.rfFreq = 12;  % MHz
ExpB.Dim1 = {'rfFreq',0.5};
ExpB.nPoints = 9;
[xB,yB,infoB] = saffron(Sys,ExpB,Opt);

common = 2:10;
ok(1) = areequal(xA(common),xB,1e-12,'abs');
ok(2) = areequal(yA(common)-yA(1),yB-yA(1),1e-6,'rel');
ok(3) = all(isfinite(yA));
ok(4) = areequal(diff(infoB.EventTimes(5,:)),rf.tp,1e-12,'abs');
