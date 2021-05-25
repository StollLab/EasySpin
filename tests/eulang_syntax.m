function ok = test()

an = rand(1,3);
R = erot(an);

eulang(R);
[a,b,c] = eulang(R);
aa = eulang(R);

ok = true;
