function ok = test()

% Test of zfield() for a simple axial case

Sys = struct('S',3/2,'D',[-1 -1 2]);
H0a = diag([3 -3 -3 3]);
H0b = zfield(Sys);

ok = areequal(H0a,H0b,1e-12,'abs');
