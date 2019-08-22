function [err,data] = test(opt,olddata)

a = planck;
b = 6.62607015e-34;
err = abs(a-b)/a>1e-9;
data = [];
