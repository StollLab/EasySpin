function [err,data] = test(opt,olddata)

a = boltzm;
b = 1.380649e-23;
err = abs(a-b)/b > 1e-10;

data = [];
