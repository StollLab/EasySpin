function [err,data] = test(opt,olddata)

a = bmagn;
b = 9.2740100783e-24;
err = abs(a-b)/b > 1e-10;
data = [];
