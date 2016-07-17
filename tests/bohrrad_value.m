function [err,data] = test(opt,olddata)

a = bohrrad;
b = 0.52917721067e-10;
err = abs(a-b)/b > 1e-10;
data = [];
