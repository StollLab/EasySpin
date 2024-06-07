function ok = test()

a = emass;
b = 9.1093837139e-31;

ok = areequal(a,b,1e-11,'rel');
