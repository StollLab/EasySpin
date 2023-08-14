function ok = test()

% Axial+D2h, collinear
Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 2 1]);
G = hamsymm(Sys);
ok(1) = strcmp(G,'D2h');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = hamsymm(Sys);
ok(2) = strcmp(G,'D2h');
