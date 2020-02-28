function ok = test()

rng(55);
v = mhz2mt;
v = mhz2mt(rand);
v = mhz2mt(rand,rand);
v = mhz2mt(rand,rand(1,6));

ok = true;
