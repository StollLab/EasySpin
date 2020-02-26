function ok = test()

N = 1000;
[x,y,z] = sphrand(N);

ok = (numel(x)==N) && (numel(y)==N) && (numel(z)==N);
