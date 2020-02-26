function ok = test()

y = rand(100,5);
a = deriv(y);

ok = areequal(size(a),size(y),0);
