function ok = test()

v = rand(3,10);
[p,t] = vec2ang(v);
a = vec2ang(v);

ok = all(a==[p;t]);
