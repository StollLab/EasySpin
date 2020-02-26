function ok = test()

% Check for a fixed bug, where saffron crashed if Sys.A_ was given
% and ExciteWidth triggered the removal of a nucleus for coreSys.

Sys.S = 1/2;
Sys.Nucs = '1H,1H';
Sys.lwEndor = 0.1;

Exp.Field = 350;
Exp.Sequence = 'MimsENDOR';

Sys.A_ = [200 20; 3 2];

Exp.tau = 444;
Exp.Range = [13 16.4];
Exp.mwFreq = 9.7;
Exp.ExciteWidth = 50;

[x,y] = saffron(Sys,Exp);

ok = true;
