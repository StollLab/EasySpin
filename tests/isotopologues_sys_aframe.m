function ok = test()

% Isotopologues with AFrame

Sys.Nucs = 'C';
Sys.n = 3;
Sys.A = [10 10 20];
Sys.AFrame = [1 2 2]*pi/3;

Iso = isotopologues(Sys,0);

ok = numel(Iso(1).AFrame)==0 && ...
     numel(Iso(2).AFrame)==3 && ...
     numel(Iso(3).AFrame)==3 && ...
     numel(Iso(4).AFrame)==3;

ok = ok && numel(Iso)==4;
