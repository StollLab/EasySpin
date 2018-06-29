% HYSCORE of a 14N nucleus (saffron)
%==========================================================================

clear, clf

Sys.Nucs = '14N';
Sys.A_ = [0.8 0.5];
Sys.Q = [-1 -1 2]*0.1;

Exp.Sequence = 'HYSCORE';
Exp.Field = 330;
Exp.tau = 0.080;
Exp.dt = 0.120;
Exp.nPoints = 256;

Opt.nKnots = 181;

saffron(Sys,Exp,Opt);
