function ok = test() %#ok<FNDEF>

% Relaxation in the manifold engine applies to any sequence: T2 during
% intervals with electron coherence, T1 during intervals with populations.
% A custom three-pulse sequence gives the same result as the predefined
% experiment, and the two-pulse echo decays as exp(-2*tau/T2) as in the
% full engine (T2 acts on electron coherences).

Sys.Nucs = '1H';
Sys.A = [2 3 9];  % MHz
Sys.T1 = 3;  % µs
Sys.T2 = 0.5;  % µs

Opt.GridSize = 7;

Exp.Field = 350;  % mT
Exp.Sequence = '3pESEEM';
Exp.tau = 0.12;  % µs
Exp.T = 0.1;  % µs
Exp.dt = 0.02;  % µs
Exp.nPoints = 64;
[~,y1,info] = saffron(Sys,Exp,Opt);
[~,y2] = saffron(Sys,rmfield(info.Exp,'predefinedExp'),Opt);
ok(1) = areequal(y1,y2,1e-12,'rel');

% two-pulse echo, single orientation, no modulation (isotropic A, high field)
S.Nucs = '1H';
S.A = 0;
S.T2 = 0.5;  % µs
p90.Flip = pi/2;
p180.Flip = pi;
E.Field = 350;  % mT
E.mwFreq = 9.8;  % GHz
E.MolFrame = [0 0 0];
E.CrystalSymmetry = 1;
E.Sequence = {p90 0.1 p180 0.1};
E.nPoints = 8;
E.Dim1 = {'d1,d2',0.05};
E.Pathways = '+-';
E.DetPhase = -pi/2;
[x,yP] = saffron(S,E,struct('Method','manifold'));
[~,yT] = saffron(S,E,struct('Method','full'));
ok(2) = areequal(real(yP),0.5*exp(-2*x/S.T2),1e-6,'rel');
ok(3) = areequal(real(yT),real(yP),1e-3,'rel');
