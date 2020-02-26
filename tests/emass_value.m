function ok = test()

a = emass;
b = 9.1093837015e-31;

ok = areequal(a,b,1e-12,'rel');
