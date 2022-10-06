function ok = test()

% Test of ham_zf() for a simple axial case

Sys = struct('S',3/2,'D',[-1 -1 2]);
H0a = diag([3 -3 -3 3]);
H0b = ham_zf(Sys);

ok(1) = areequal(H0a,H0b,1e-12,'abs');

clear Sys
Sys = struct('S',3/2,'D_',[3 0]);
H0a = diag([3 -3 -3 3]);
H0b = ham_zf(Sys);

ok(2) = areequal(H0a,H0b,1e-12,'abs');