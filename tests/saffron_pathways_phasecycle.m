function ok = test() %#ok<FNDEF>

% Pathway selection with Exp.Pathways is equivalent to an ideal phase
% cycle (Exp.PhaseCycle), in both engines: three-pulse ESEEM, single
% orientation, ideal pulses.

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz

p90.Flip = pi/2;

Exp.Field = 350;  % mT
Exp.mwFreq = 9.8;  % GHz
Exp.MolFrame = [0 0.4 0.5];
Exp.CrystalSymmetry = 1;
Exp.Sequence = {p90 0.12 p90 0.3 p90 0.12};
Exp.DetPhase = -pi/2;

sel = [0 1; pi/2 1i; pi -1; 3*pi/2 -1i];
odd = [0 1; pi -1];

E1 = Exp; E1.Pathways = '+0-';
E2 = Exp; E2.PhaseCycle = {sel [] odd};

tol = 1e-6;
methods = {'subspace','fullspace'};
for m = 1:2
  Opt.Method = methods{m};
  [~,y1] = saffron(Sys,E1,Opt);
  [~,y2] = saffron(Sys,E2,Opt);
  ok(m) = areequal(y1,y2,tol,'rel');
end
