function ok = test()

% Test Two Spin: Hyperfine Anisotropy

Sys = struct('S',[1/2 1/2],'g',[2 2],'Nucs','1H','A',[0 0 0 1 1 2],'ee',0);
G = symm(Sys);
ok = strcmp(G,'Dinfh');

Sys.ZB11.l = 0;
Sys.ZB11.vals = 1;
G = symm(Sys);
ok = strcmp(G,'Dinfh');
