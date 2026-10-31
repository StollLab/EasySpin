function ok = test()

% Relaxation in thyme mode
% (1) matrix T2 refers to product states: with isotropic hyperfine coupling,
%     the pulses only create the allowed electron coherences between product
%     states 1-3 and 2-4 (aa-ba, ab-bb), so T2 for these two gives the same
%     echo as a scalar T2. With A > 2*nuI, the eigenstates sorted by energy
%     are in a different order than the product states.
% (2) short T1 does not lead to numerical divergence
% (3) scalar T2 acts only on electron coherences: in three-pulse ESEEM with
%     anisotropic hyperfine coupling, it gives the same result as T2 for
%     all electron coherences (aa,ab - ba,bb), and nuclear coherences
%     during T do not decay

Sys.S = 1/2;
Sys.g = 2;
Sys.Nucs = '1H';
Sys.A = 40;  % MHz

Exp.Field = 350;
Exp.mwFreq = 9.8;
P90.Type = 'rectangular';
P90.tp = 0.008;
P90.Flip = pi/2;
P180 = P90;
P180.tp = 0.016;
P180.Flip = pi;
Exp.Sequence = {P90 0.2 P180 0.2};
Exp.DetWindow = [-0.004 0.004];

Opt.GridSize = 7;

% (1)
[~,y0] = saffron(Sys,Exp,Opt);
SysR = Sys;
SysR.T2 = 0.3;
[~,yScalar] = saffron(SysR,Exp,Opt);
SysR.T2 = zeros(4);
SysR.T2(1,3) = 0.3;
SysR.T2(2,4) = 0.3;
[~,yMatrix] = saffron(SysR,Exp,Opt);
ok(1) = areequal(yMatrix,yScalar,1e-8,'rel');
ok(2) = abs(sum(yScalar)) < 0.5*abs(sum(y0));  % echo is damped

% (2) S = 1/2 with 14N (6 states), T1 much shorter than the sequence
Sys.Nucs = '14N';
Sys.A = [10 30];
Exp.Sequence = {P180 1.5 P90 0.2 P180 0.2};
[~,y0] = saffron(Sys,Exp,Opt);
Sys.T1 = 0.03;
[~,y] = saffron(Sys,Exp,Opt);
ok(3) = all(isfinite(y)) && max(abs(y)) < 2*max(abs(y0));

% (3)
Sys = rmfield(Sys,'T1');
Sys.Nucs = '1H';
Sys.A = [2 9];
Exp = rmfield(Exp,'DetWindow');
Exp.Sequence = '3pESEEM';
Exp.tau = 0.12;
Exp.T = 0.1;
Exp.dt = 0.1;
Exp.nPoints = 10;
Opt.SimulationMode = 'thyme';
SysR = Sys;
SysR.T2 = 0.5;
[~,yScalar] = saffron(SysR,Exp,Opt);
SysR.T2 = zeros(4);
SysR.T2([1 2],[3 4]) = 0.5;
[~,yMatrix] = saffron(SysR,Exp,Opt);
ok(4) = areequal(yScalar,yMatrix,1e-8,'rel');
