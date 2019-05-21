function [err,data] = test(opt,olddata)

a = rydberg;
b = 10973731.568160;
err = abs(a-b)/b > 1e-12;

data = [];
