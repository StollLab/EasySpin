function [err,data] = test(opt,olddata)

a = planck;
b = 6.626070040e-34;
err = abs(a-b)/a>1e-10;
data = [];
