function [err,data] = test(opt,olddata)

a = molgas;
b = 8.3144598;
err = ~areequal(a,b);
data = [];
