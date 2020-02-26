function ok = test()

a = clebschgordan(0,0,0,0,0,0);
b = 1;
ok = areequal(a,b,1e-10,'rel');

