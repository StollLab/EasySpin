% HYSCORE of a 14N nucleus (saffron)
%==========================================================================

clear, clf

Sys.Nucs = '14N';
Sys.A_ = [0.8 0.5]; % MHz
Sys.Q = [-1 -1 2]*0.1; % MHz

Exp.Sequence = 'HYSCORE';
Exp.Field = 330; % mT
Exp.tau = 0.080; % mus
Exp.dt = 0.120; % mus
Exp.nPoints = 256;

Opt.nKnots = 181;

saffron(Sys,Exp,Opt);
