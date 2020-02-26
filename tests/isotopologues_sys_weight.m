function ok = test()

% Test whether isotopologue() included given Sys.weight

w = 0.8;

Sys.g = 2;
Sys.Nucs = '(1,2)H';
Sys.A = 10;
Sys.Abund = [w 1-w];
Sys.weight = 0.1;

Iso = isotopologues(Sys);

ok = Iso(1).weight==w*Sys.weight && Iso(2).weight==(1-w)*Sys.weight;
