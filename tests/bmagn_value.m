function [err,data] = test(opt,olddata)

a = bmagn;
b = 9.274009994e-24;
err = abs(a-b)/b > 1e-10;
data = [];
