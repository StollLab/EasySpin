function [err,data] = test(opt,olddata)

% Assert that the predefined and user-defined HYSCORE experiments give the same result

Sys.Nucs = '14N';
Sys.A = [5 6 7];
Sys.Q = 2;

Exp.CrystalOrientation = [0 rand*pi 0];

Exp.Field = 350;

tau = 0.100;
dt = 0.030; % time increment

Exp.Flip = [1 1 2 1];
Exp.Inc = [0 1 2 0];
Exp.t = [tau 0 0 tau];
Exp.dt = dt;
Opt.TimeDomain = true;

[x,y] = saffron(Sys,Exp,Opt);

Exp2.CrystalOrientation = Exp.CrystalOrientation;
Exp2.Field = Exp.Field;
Exp2.Sequence = 'HYSCORE';
Exp2.dt = dt;
Exp2.tau = tau;
Opt.TimeDomain = true;

[x,y2] = saffron(Sys,Exp2,Opt);
y2 = -y2;

err = max(abs(y(:)-y2(:)))>1e-10;

data = [];
