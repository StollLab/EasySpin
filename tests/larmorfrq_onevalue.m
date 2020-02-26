function ok = test()

Nuc = '1H';
B = 100;
nu1 = larmorfrq(Nuc,B);
nu2 = abs(nucgval(Nuc))*B*nmagn/planck/1e9;
ok = areequal(nu1,nu2,1e-10,'rel');
