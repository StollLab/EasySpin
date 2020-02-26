function ok = test()

N = 400;
v = sphrand(N);
ok = all(size(v)==[3 N]);

