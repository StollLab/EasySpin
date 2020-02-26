function ok = test()

a = faraday;
b = avogadro*echarge;
ok = areequal(a,b,1e-10,'rel');
