function ok = test()

Sys1.g = 2;
Sys1.Nucs = 'Cl';
Sys1.A = [2 5];

Sys2.g = 2;
Sys2.Nucs = 'Cl,Cl';
Sys2.A = [2 5; 1 3];

Exp.Sequence = '2pESEEM';
Exp.Field = 350;
Exp.dt = 0.010;
Exp.nPoints = 200;

Opt.Verbosity = 0;
Opt.separate = 'components';

[t,V] = saffron(Sys1,Exp,Opt);
ok(1) = size(V,1)==2;  % 1 component with 2 isotopologues

[t,V] = saffron(Sys2,Exp,Opt);
ok(2) = size(V,1)==4;  % 1 component with 4 isotopologues
