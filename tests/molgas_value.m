function ok = test()

a = molgas;
b = avogadro*boltzm;
ok = areequal(a,b,1e-10,'rel');
