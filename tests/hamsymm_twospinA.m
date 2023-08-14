function ok = test()

% Test Two Spin: Hyperfine Anisotropy

Sys = struct('S',[1/2 1/2],'g',[2 2],'Nucs','1H','A',[1 1 1 2 2 3],'ee',6);
G = hamsymm(Sys);
ok(1) = strcmp(G,'Dinfh');

Sys.AFrame = [0 0 0 0 0 0];
G = hamsymm(Sys);
ok(2) = strcmp(G,'Dinfh');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = hamsymm(Sys);
ok(3) = strcmp(G,'Dinfh');
