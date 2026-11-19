function ok = test()  %#ok<FNDEF>

% Exp.PhaseCycle works in both engines. With finite pulses, the pathways
% engine weights each coherence pathway with the cycle, so the echo-selecting
% cycle gives the same result as selecting the echo pathway with
% Exp.Pathways, and the full-space engine agrees.

Sys.S = 1/2;
Sys.g = gfree;
Sys.Nucs = '1H';
Sys.A = [2 6];  % MHz

p90 = struct('tp',0.016,'Flip',pi/2);
p180 = struct('tp',0.032,'Flip',pi);
tau = 0.1;  % µs

Exp.mwFreq = 9.5;  % GHz
Exp.Field = unitconvert(Exp.mwFreq,'GHz->mT',gfree);
Exp.MolFrame = [0 0 0];
Exp.SampleFrame = [0 pi/4 0];
Exp.Sequence = {p90 tau p180 tau};
Exp.Dim1 = {'d1,d2' 0.02};
Exp.nPoints = 3;

ph = (0:3).'*pi/2;
E1 = Exp; E1.PhaseCycle = [ph 0*ph -ph];
E2 = Exp; E2.Pathways = '+-';

[~,y1,info] = saffron(Sys,E1);
[~,y2] = saffron(Sys,E2);
ok(1) = strcmp(info.Method,'subspace');
ok(2) = areequal(y1,y2,1e-10,'rel');

[~,y3] = saffron(Sys,E1,struct('Method','fullspace'));
ok(3) = areequal(y3,y1,1e-3,'rel');
