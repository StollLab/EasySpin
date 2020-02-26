function ok = test()

a = rand(1,3);
R1 = erot(a);
b = eulang(R1);
ok = areequal(a,b,1e-10,'abs');
