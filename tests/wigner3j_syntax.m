function ok = test()

% Compare outputs for all calling syntaxes

J1 = 4;
J2 = 3;
J3 = 2;
M1 = 1;
M2 = -2;
M3 = 1;

a0 = -sqrt(7/5)/6;

a1 = wigner3j(J1,J2,J3,M1,M2,M3);
a2 = wigner3j([J1 J2 J3],[M1 M2 M3]);
a3 = wigner3j([J1 M1],[J2 M2],[J3 M3]);
a4 = wigner3j([J1 J2 J3; M1 M2 M3]);

thr = 1e-10;
ok(1) = areequal(a1,a0,thr,'abs');
ok(2) = areequal(a2,a0,thr,'abs');
ok(3) = areequal(a3,a0,thr,'abs');
ok(4) = areequal(a4,a0,thr,'abs');
