function ok = test() %#ok<FNDEF>

% Pathway selection with Exp.Pathways is equivalent to an ideal phase
% cycle (Exp.PhaseCycle), in both engines: three-pulse ESEEM, single
% orientation, ideal pulses.

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

% Explicit phase cycles per pulse
pc1 = [0 1; pi/2 1i; pi -1; 3*pi/2 -1i];
pc2 = [];
pc3 = [0 1; pi -1];

Exp_string = Exp;
Exp_string.Pathways = '+0-';
Exp_cell = Exp;
Exp_cell.PhaseCycle = {pc1 pc2 pc3};

tol = 1e-6;
methods = {'subspace','fullspace'};
for m = 1:2
  Opt.Method = methods{m};
  [~,y1] = saffron(Sys,Exp_string,Opt);
  [~,y2] = saffron(Sys,Exp_cell,Opt);
  ok(m) = areequal(y1,y2,tol,'rel');
end
