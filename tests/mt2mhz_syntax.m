function ok = test()

rng(823);
v = mt2mhz;
v = mt2mhz(rand);
v = mt2mhz(rand,rand);
v = mt2mhz(rand,rand(1,6));

ok = true;
