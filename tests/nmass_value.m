function [err,data] = test(opt,olddata)

a = nmass;
b = 1.674927471e-27;
err = abs(a-b)/b > 1e-12;

data = [];
