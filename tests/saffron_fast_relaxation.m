function ok = test()

% Fast mode: relaxation decays multiplied onto predefined sequences
% (T2 over 2*tau, T1 over the population periods)

Sys.Nucs = '1H';
Sys.A = [2 9];

Exp.Field = 350;
Exp.tau = 0.1;
Exp.T = 0.08;
Exp.dt = 0.01;
Exp.nPoints = 50;

Opt.GridSize = 7;

T1 = 2;
T2 = 0.5;
t = (0:Exp.nPoints-1)*Exp.dt;
tol = 1e-10;

% Two-pulse ESEEM, T2 over 2*(tau+t)
Exp.Sequence = '2pESEEM';
[~,y0] = saffron(Sys,Exp,Opt);
SysR = Sys; SysR.T2 = T2;
[~,y] = saffron(SysR,Exp,Opt);
ok(1) = areequal(y,y0.*exp(-2*(Exp.tau+t)/T2),tol,'rel');

% Three-pulse ESEEM with only T2: constant factor
Exp.Sequence = '3pESEEM';
[~,y0] = saffron(Sys,Exp,Opt);
[~,y] = saffron(SysR,Exp,Opt);
ok(2) = areequal(y,y0*exp(-2*Exp.tau/T2),tol,'rel');

% Three-pulse ESEEM, T1 over T+t
SysR.T1 = T1;
[~,y] = saffron(SysR,Exp,Opt);
ok(3) = areequal(y,y0*exp(-2*Exp.tau/T2).*exp(-(Exp.T+t)/T1),tol,'rel');

% Four-pulse ESEEM, T1 over 2*(T+t)
Exp.Sequence = '4pESEEM';
[~,y0] = saffron(Sys,Exp,Opt);
[~,y] = saffron(SysR,Exp,Opt);
ok(4) = areequal(y,y0*exp(-2*Exp.tau/T2).*exp(-2*(Exp.T+t)/T1),tol,'rel');

% T1 = T2 = 0: no relaxation
SysR.T1 = 0;
SysR.T2 = 0;
[~,y] = saffron(SysR,Exp,Opt);
ok(5) = areequal(y,y0,tol,'rel');
