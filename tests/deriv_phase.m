function ok = test()

y = 1:100;
a = deriv(y);

ok = all(a>0);
