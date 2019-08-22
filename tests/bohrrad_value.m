function [err,data] = test(opt,olddata)

a = bohrrad;
b = 0.529177210903e-10;
err = abs(a-b)/b > 1e-10;
data = [];
