function ok = test()

% D2h+D2h, general tilt

Sys = struct('S',1/2,'Nucs','1H','g',[2 2.5 3],'A',[3 2 1],'AFrame',rand(1,3)+pi/30);
G = symm(Sys);
ok(1) = strcmp(G,'Ci');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = symm(Sys);
ok(2) = strcmp(G,'Ci');
