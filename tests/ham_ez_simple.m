function ok = test()

Sys = struct('S',1/2,'g',[2 2.1 2.2]);

[mux0,muy0,muz0] = ham_ez(Sys);

f = -Sys.g*bmagn/planck/1e9;  % MHz/mT
mux1 = f(1)*sop(Sys,'x');
muy1 = f(2)*sop(Sys,'y');
muz1 = f(3)*sop(Sys,'z');

thr = 1e-10;
ok(1) = areequal(mux0,mux1,thr,'abs');
ok(2) = areequal(muy0,muy1,thr,'abs');
ok(3) = areequal(muz0,muz1,thr,'abs');
