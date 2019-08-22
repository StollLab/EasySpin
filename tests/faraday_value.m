function [err,data] = test(opt,olddata)

a = faraday;
b = avogadro*echarge;
err = abs(a-b)/b > 1e-10;
data = [];
