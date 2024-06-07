function ok = test()

a = pmass;
b = 1.67262192595e-27;
ok = areequal(a,b,1e-12,'rel');
