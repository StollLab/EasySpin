function [err,data] = test(opt,olddata)

a = emass;
b = 9.10938356e-31;
err = abs(a-b)/b > 1e-12;
data = [];