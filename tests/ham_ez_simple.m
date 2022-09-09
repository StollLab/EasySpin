function ok = test()

Sys = struct('S',1/2,'g',[2 2.1 2.2]);
[Gx0,Gy0,Gz0] = ham_ez(Sys);
f = -Sys.g*bmagn/planck/1e9;  % MHz/mT

mux1 = f(1)*sop(Sys,'x');
muy1 = f(2)*sop(Sys,'y');
muz1 = f(3)*sop(Sys,'z');

A = [Gx0 Gy0 Gz0];
B = [mux1 muy1 muz1];

ok = areequal(A,B,1e-10,'abs');
