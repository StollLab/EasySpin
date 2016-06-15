function [err,data] = test(opt,olddata)

% Test 1:
%======================================================
tri = sphtri('D2h',2);
err = any(tri~=[1 2 3]);
data = [];
