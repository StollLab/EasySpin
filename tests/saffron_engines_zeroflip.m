function ok = test() %#ok<FNDEF>

% A finite microwave pulse with Flip = 0 is free evolution over its length,
% in both engines. The difference to the signal with the pulse on removes
% all signal that does not depend on the pulse (as with a pulse switched on
% and off on a spectrometer); both engines agree on it.

Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = [3 8];  % MHz

Exp.Field = 1200;  % mT
Exp.mwFreq = mt2mhz(Exp.Field,Sys.g)/1e3;  % GHz, on resonance
Exp.SampleFrame = [0 -0.7 0];
pInv = struct('Flip',pi,'tp',0.032);  % µs
p90 = struct('Flip',pi/2,'tp',0.016);
p180 = struct('Flip',pi,'tp',0.032);
Exp.Sequence = {pInv 0.4 p90 0.3 p180 0.3};

Eoff = Exp;
Eoff.Sequence{1}.Flip = 0;

methods = {'subspace','fullspace'};
for m = 1:2
  Opt.Method = methods{m};
  [~,yOn(m)] = saffron(Sys,Exp,Opt); %#ok<AGROW>
  [~,yOff(m)] = saffron(Sys,Eoff,Opt); %#ok<AGROW>
end
yDiff = yOn - yOff;

% without the inversion pulse: the two-pulse echo after a delay of
% 0.4 µs + tp, which does not change the echo of an isolated spin system
E2p = Exp;
E2p.Sequence = {p90 0.3 p180 0.3};
[~,y2p] = saffron(Sys,E2p,struct('Method','subspace'));

ok(1) = areequal(yOff(1),yOff(2),1e-4,'rel');
ok(2) = areequal(yDiff(1),yDiff(2),1e-4,'rel') && abs(yDiff(1))>0.1*abs(yOff(1));
ok(3) = areequal(yOff(1),y2p,1e-10,'rel');
