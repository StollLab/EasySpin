function ok = test()

% Two axial, collinear tensors

Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 3 1]);
G = symm(Sys);
ok(1) = strcmp(G,'Dinfh');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = symm(Sys);
ok(2) = strcmp(G,'Dinfh');
