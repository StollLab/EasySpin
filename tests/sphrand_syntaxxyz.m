function [err,data] = test(opt,olddata)

% Test 2: Calling syntax: phi, theta
%=======================================================
N = 1000;
[x,y,z] = sphrand(N);
err = (numel(x)~=N) || (numel(y)~=N) || (numel(z)~=N);

data = [];
