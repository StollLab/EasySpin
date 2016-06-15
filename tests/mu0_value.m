function [err,data] = test(opt,olddata)

a = mu0;
b = 4*pi*1e-7;

err = abs(a-b)/a>1e-10;

data = [];
