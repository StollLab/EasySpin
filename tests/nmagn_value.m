function [err,data] = test(opt,olddata)

a = nmagn;
b = 5.050783699e-27;
err = abs(a-b)/b > 1e-10;
data = [];
