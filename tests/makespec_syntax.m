function ok = test()

s1 = makespec([0 1],1000,0.5,1);
[x,s2] = makespec([0 1],1000,0.5,1);

ok = all(abs(s1-s2))<1e-10;
