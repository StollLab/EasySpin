function [err,data] = test(opt,olddata)

a = rydberg;
b = 10973731.568508;
err = abs(a-b)/b > 1e-12;

data = [];
