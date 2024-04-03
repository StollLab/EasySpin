function ok = test()

% 2x Axial, z axis tilt 90 deg

Sys = struct('S',1/2,'Nucs','1H','g',[2 2 3],'A',[3 3 1],'AFrame',[0 -pi/2 0]);
G = hamsymm(Sys);
ok(1) = strcmp(G,'D2h');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = hamsymm(Sys);
ok(2) = strcmp(G,'D2h');
