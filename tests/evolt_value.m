function [err,data] = test(opt,olddata)

a = evolt;
b = echarge;
err = abs(a-b)/b > 1e-10;
data = [];
