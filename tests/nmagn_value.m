function [err,data] = test(opt,olddata)

a = nmagn;
b = 5.0507837461e-27;
err = abs(a-b)/b > 1e-11;
data = [];
