function [err,data] = test(opt,olddata)

% Direct value check

a = gfree;
b = 2.00231930436256;
err = abs(a-b)/b > 1e-14;
data = [];
