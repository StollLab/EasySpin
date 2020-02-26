function ok = test()

x_mT = 1;
g = 2;
v1 = mt2mhz(x_mT,g);
v2 = x_mT*1e-3*g*bmagn/planck/1e6;

ok = areequal(v1,v2,1e-10,'rel');
