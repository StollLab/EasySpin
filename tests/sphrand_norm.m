function [err,data] = test(opt,olddata)

% test 4: Normalized vectors?
%=======================================================
N = 10000;
v = sphrand(N);
norms = sum(v.^2);
err = max(abs(norms-1))>1e-10;
data = [];
