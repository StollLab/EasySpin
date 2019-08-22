function [err,data] = test(opt,olddata)

a = avogadro;
b = 6.02214076e23;
err = abs(a-b)/b > 1e-10;

data = [];
