function ok = test()

% Isotopologues with AFrame

Sys.Nucs = 'Cu';
Sys.n = [3];
Sys.A = [10 10 20];
Sys.Q = [10];
Sys.QFrame = [1 2 2]*pi/3;

Iso = isotopologues(Sys);

ok = numel(Iso(1).QFrame)==3 && ...
     numel(Iso(2).QFrame)==6 && ...
     numel(Iso(3).QFrame)==6 && ...
     numel(Iso(4).QFrame)==3;
