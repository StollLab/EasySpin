function [err,data] = test(opt,olddata)

a = pmass;
b = 1.67262192369e-27;
err = abs(a-b)/b > 1e-12;
data = [];
