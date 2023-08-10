function ok = test()

% Isotropic, tilted frame

Sys = struct('S',1/2,'g',[2 2 2],'gFrame',rand(1,3));
G = hamsymm(Sys);
ok(1) = strcmp(G,'O3');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = hamsymm(Sys);
ok(2) = strcmp(G,'O3');
