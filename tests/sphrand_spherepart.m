function ok = test()

N = 1000;
[x,y,z] = sphrand(N,4);
ok(1) =  min(z)>=0;

[x,y,z] = sphrand(N,2);
ok(2) = min(y)>=0 && min(z)>=0;

[x,y,z] = sphrand(N,1);
ok(3) = min(y)>=0 && min(z)>=0 && min(x)>=0;
