function [err,data] = test(opt,olddata)

% Test 2: Calling syntax: phi, theta
%=======================================================
N = 1000;
[p,t] = sphrand(N);
err = (numel(p)~=N) | (numel(t)~=N);
data = [];
