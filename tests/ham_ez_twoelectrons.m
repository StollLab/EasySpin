function ok = test()

% Two electrons

Sys = struct('S',[1/2 1],'g',[2 2.1 2.2; 3 4 5],'ee',[1 2 3]);
[mux0,muy0,muz0] = ham_ez(Sys,2);
f = -Sys.g(2,:)*bmagn/planck/1e9;

mux1 = f(1)*sop(Sys,'ex');
muy1 = f(2)*sop(Sys,'ey');
muz1 = f(3)*sop(Sys,'ez');

A = [mux0 muy0 muz0];
B = [mux1 muy1 muz1];

ok = areequal(A,B,1e-10,'abs');
