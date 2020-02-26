function ok = test()

Nuc = '14N';
B = 100:200;
nu1 = larmorfrq(Nuc,B);
nu2 = abs(nucgval(Nuc))*B*nmagn/planck/1e9;
ok = areequal(nu1,nu2.',1e-10,'rel');
