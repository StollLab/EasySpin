function [err,data] = test(opt,olddata)

Sys.Nucs = '14N';
Sys.A = [-1 -1 2]*0.5+0.8;
Sys.Q = 1;

% User-defined HYSCORE experiment
tau = 0.080;  % us
Exp.Flip = [1 1 2 1];
Exp.t = [tau 0 0 tau];
Exp.Inc = [0 1 2 0];
Exp.dt = 0.025;

Exp.Field = 330;
Exp.tau = 0.001;
Exp.dt = 0.1;
Exp.nPoints = 200;

Opt.Verbosity = 0;
y = saffron(Sys,Exp,Opt);

data = [];
err = false;
