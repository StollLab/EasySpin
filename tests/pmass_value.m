function [err,data] = test(opt,olddata)

a = pmass;
b = 1.672621898e-27;
err = abs(a-b)/b > 1e-12;
data = [];
