function [err,data] = test(opt,olddata)

% Test 1: Calling syntax: vectors
%=======================================================
N = 400;
v = sphrand(N);
err = any(size(v)~=[3 N]);
data = [];
