function [err,data] = test(opt,olddata)

a = faraday;
b = 96485.33289;
err = abs(a-b)/b > 1e-10;
data = [];
