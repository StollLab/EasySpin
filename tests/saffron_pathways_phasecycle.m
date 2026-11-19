function ok = test() %#ok<FNDEF>

% Pathway selection with Exp.Pathways is equivalent to an ideal phase
% cycle (Exp.PhaseCycle), in both engines: three-pulse ESEEM, single
% orientation, ideal pulses. The nested cycle (first pulse in four steps of
% 90°, third pulse in two steps of 180°) selects the stimulated echo. So
% does the nested cycle extended by CYCLOPS (all pulses and the receiver
% stepped together by 90°), which is not a nested cycle.

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

p90.Flip = pi/2;  % rad

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0.4 0.5];
Exp.CrystalSymmetry = 1;
tau = 0.12;  % µs
T = 0.3;  % µs
Exp.Sequence = {p90 tau p90 T p90 tau};
Exp.DetPhase = -pi/2;  % rad

% Phase cycles, rows [phase1 phase2 phase3 receiverPhase]
[ph1,ph3] = ndgrid((0:3)*pi/2,[0 pi]);
nested = [ph1(:) 0*ph1(:) ph3(:) -ph1(:)-ph3(:)];
[iStep,iCyclops] = ndgrid(1:8,0:3);
cyclops = nested(iStep(:),:) + iCyclops(:)*pi/2*[1 1 1 1];

Exp_pathways = Exp;
Exp_pathways.Pathways = '+0-';
Exp_nested = Exp;
Exp_nested.PhaseCycle = nested;
Exp_cyclops = Exp;
Exp_cyclops.PhaseCycle = cyclops;

tol = 1e-6;
methods = {'subspace','fullspace'};
for m = 1:2
  Opt.Method = methods{m};
  [~,y1] = saffron(Sys,Exp_pathways,Opt);
  [~,y2] = saffron(Sys,Exp_nested,Opt);
  [~,y3] = saffron(Sys,Exp_cyclops,Opt);
  ok(2*m-1) = areequal(y1,y2,tol,'rel'); %#ok<AGROW>
  ok(2*m) = areequal(y1,y3,tol,'rel'); %#ok<AGROW>
end
