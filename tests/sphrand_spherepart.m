function [err,data] = test(opt,olddata)

err = 0;

N = 1000;
[x,y,z] = sphrand(N,4);
err = err || min(z)<0;

[x,y,z] = sphrand(N,2);
err = err || min(y)<0 || min(z)<0;

[x,y,z] = sphrand(N,1);
err = err || min(y)<0 || min(z)<0 || min(x)<0;

data = [];
