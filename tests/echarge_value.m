function [err,data] = test(opt,olddata)

a = echarge;
b = 1.6021766208e-19;
err = abs(a-b)/b > 1e-10;
data = [];
