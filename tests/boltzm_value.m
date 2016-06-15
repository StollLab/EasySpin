function [err,data] = test(opt,olddata)

a = boltzm;
b = 1.38064852e-23;
err = abs(a-b)/b > 1e-10;

data = [];
