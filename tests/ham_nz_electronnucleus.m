function ok = test()

Nuc = '14N';
Sys = struct('S',1/2,'Nucs',Nuc,'g',[2 2.1 2.2],'A',[1 2 3]);
[mux0,muy0,muz0] = ham_nz(Sys,1);
pre = +nucgval(Nuc)*nmagn/planck/1e9;

mux1 = pre*sop(Sys,'ex');
muy1 = pre*sop(Sys,'ey');
muz1 = pre*sop(Sys,'ez');

thr = 1e-10;
ok(1) = areequal(mux0,mux1,thr,'abs');
ok(2) = areequal(muy0,muy1,thr,'abs');
ok(3) = areequal(muz0,muz1,thr,'abs');
