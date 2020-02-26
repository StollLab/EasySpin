function ok = test()

N = 1000;
[p,t] = sphrand(N);
ok = numel(p)==N && numel(t)==N;
