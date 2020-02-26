function ok = test()

a = clebschgordan([300 0],[300 1],[300 1]);
b = 0.0247298820065665620596214763779;

ok = areequal(a,b,1e-10,'rel');
