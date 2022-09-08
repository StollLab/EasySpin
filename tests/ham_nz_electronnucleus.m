function ok = test()

Nuc = '14N';
Sys = struct('S',1/2,'Nucs',Nuc,'g',[2 2.1 2.2],'A',[1 2 3]);
[Gx0,Gy0,Gz0] = ham_nz(Sys,1);
pre = -nucgval(Nuc)*nmagn/planck/1e9;

Gx1 = pre*sop(Sys,'ex');
Gy1 = pre*sop(Sys,'ey');
Gz1 = pre*sop(Sys,'ez');

thr = 1e-10;
ok(1) = areequal(Gx0,Gx1,thr,'abs');
ok(2) = areequal(Gy0,Gy1,thr,'abs');
ok(3) = areequal(Gz0,Gz1,thr,'abs');
