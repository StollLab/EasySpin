function ok = test()

% Axial+D2h, general tilt

Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 2 1],'gFrame',rand(1,3)+pi/10);
G = hamsymm(Sys);
ok(1) = strcmp(G,'Ci');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = hamsymm(Sys);
ok(2) = strcmp(G,'Ci');
