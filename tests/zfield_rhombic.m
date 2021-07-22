function ok = test()

% Test of zfield() for a simple rhombic case

Sys = struct('S',1,'D',[3 1]);
H0a = [1 0 1; 0 -2 0; 1 0 1];
H0b = zfield(Sys);
ok(1) = areequal(H0a,H0b,1e-12,'abs');

Sys = struct('S',1,'D_',[3 1/3]);
H0a = [1 0 1; 0 -2 0; 1 0 1];
H0b = zfield(Sys);
ok(2) = areequal(H0a,H0b,1e-12,'abs');