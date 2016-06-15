function [err,data] = test(opt,olddata)

a = amu;
b = 1.660539040e-27;

err = abs(a-b)/b > 1e-12;

data = [];
