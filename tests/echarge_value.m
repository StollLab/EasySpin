function [err,data] = test(opt,olddata)

a = echarge;
b = 1.602176634e-19;
err = abs(a-b)/b > 1e-10;
data = [];
