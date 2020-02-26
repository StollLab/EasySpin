function ok = test()

val = avogadro;
ref = 6.02214076e23;
ok = areequal(val,ref,1e-10,'rel');
