function ok = test()

Sys = struct('S',1/2,'g',[2 2.1 2.2]);
[Gx0,Gy0,Gz0] = ham_ez(Sys);
f = Sys.g*bmagn/planck/1e9;  % MHz/mT

Gx1 = f(1)*sop(Sys,'x');
Gy1 = f(2)*sop(Sys,'y');
Gz1 = f(3)*sop(Sys,'z');

A = [Gx0 Gy0 Gz0];
B = [Gx1 Gy1 Gz1];

ok = areequal(A,B,1e-10,'abs');
