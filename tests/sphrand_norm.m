function ok = test()

N = 10000;
v = sphrand(N);
norms = sum(v.^2);
ok = max(abs(norms-1))<1e-10;

